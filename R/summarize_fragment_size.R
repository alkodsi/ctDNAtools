#' Summarizes fragment size in defined genomic regions
#'
#' @param bam the input bam file
#' @param regions data frame containing the genomic regions. Must have the columns chr, start and end.
#' @param tag  the RG tag if the bam has more than one sample.
#' @param summary_functions a named list containing the R functions used for summarization, e.g. mean, sd.
#' @param ... Other parameters passed to \code{\link{get_fragment_size}}
#' @return a data frame with the first column having the regions in the format of chr:start-end, and other columns correspond to summary_functions.

#' @details Fragment size for reads that are paired (optionally properly paired), whose both mates are mapped,
#' not secondary or supplementary alignment, not duplicates, passed quality control, and satisfy mapq threshold will be used for summarization.
#' The reads that overlap the specified regions will be summarized by the specified summary_functions. Overlaps consider fragments to span the
#' left most to the right most coordinate from either the read or the mate. Minimum and maximum bounds of the fragment size will be applied before summarization.

#' @seealso \code{\link{get_fragment_size}} \code{\link{bin_fragment_size}} \code{\link{analyze_fragmentation}}

#' @export
#' @importFrom rlang .data
#' @importFrom magrittr %>%
#' @importFrom stats median
#' @examples
#' data("targets", package = "ctDNAtools")
#' bamT1 <- system.file("extdata", "T1.bam", package = "ctDNAtools")
#'
#' ## binning the target in arbitrary way
#' ## Note that regions don't need to be bins,
#' ## they can be any regions in the genome
#' regions <- data.frame(
#'   chr = targets$chr,
#'   start = seq(from = targets$start - 200, to = targets$end + 200, by = 30),
#'   stringsAsFactors = FALSE
#' )
#' regions$end <- regions$start + 50
#'
#' ## basic usage
#' sfs <- summarize_fragment_size(bam = bamT1, regions = regions)
#'
#' ## different summary functions
#' sfs <- summarize_fragment_size(
#'   bam = bamT1, regions = regions,
#'   summary_functions = list(
#'     Var = var, SD = sd,
#'     meanSD = function(x) mean(x) / sd(x)
#'   )
#' )
summarize_fragment_size <- function(bam, regions, tag = "",
                                    summary_functions = list(Mean = mean, Median = median), ...) {
  assertthat::assert_that(
    assertthat::has_name(regions, c("chr", "start", "end")),
    is.data.frame(regions), assertthat::not_empty(regions)
  )

  assertthat::assert_that(!missing(bam), is.character(bam), length(bam) == 1, file.exists(bam))

  assertthat::assert_that(is.character(tag), length(tag) == 1)

  assertthat::assert_that(is.list(summary_functions), !is.null(names(summary_functions)),
    msg = "summary_functions must be a named list"
  )

  assertthat::assert_that(all(purrr::map_lgl(summary_functions, is.function)),
    msg = "all elements of summary_functions must be functions"
  )

  assertthat::assert_that(all(purrr::map_dbl(purrr::invoke_map(summary_functions, x = c(1:5)), length) == 1),
    msg = "Functions in summary_functions must produce output of length 1"
  )

  sm <- get_bam_SM(bam = bam, tag = tag)

  ellipsis::check_dots_used()

  reads <- get_fragment_size(bam = bam, tag = tag, targets = regions, ...)

  reads_gr <- GenomicRanges::GRanges(reads$chr, IRanges::IRanges(reads$start, reads$end))

  regions$ID <- paste0(regions$chr, ":", regions$start, "-", regions$end)

  regions_gr <- GenomicRanges::GRanges(regions$chr, IRanges::IRanges(regions$start, regions$end))

  overlaps <- as.data.frame(GenomicRanges::findOverlaps(reads_gr, regions_gr))

  fragments <- reads %>%
    dplyr::mutate(Region = NA) %>%
    dplyr::select(Fragment = .data$ID, .data$size, .data$Region)

  fragments[overlaps[, 1], "Region"] <- regions[overlaps[, 2], "ID"]

  summary <- fragments %>%
    tidyr::replace_na(list(Region = "offTargets")) %>%
    dplyr::group_by(.data$Region) %>%
    dplyr::summarize_if(is.numeric, summary_functions) %>%
    as.data.frame()

  colnames(summary)[-1] <- paste(sm, names(summary_functions), sep = "_")

  return(summary)
}
