#' Provides fragment ends analysis
#'
#' Calculates the number of fragment ends and the Windowed Protection Score (WPS) in genomic tiles within targets
#'
#' @param bam the input bam file
#' @param targets The targets to restrict the windows within. Must have the columns chr, start and end.
#'  In case of whole-genome, specify full chromosomes targets.
#' @param tag  the RG tag if the bam has more than one sample.
#' @param window_size The window (bin) size to use within the targets
#' @param step_size The step size to use in case of overlapping bins.
#' @param min_size Restrict fragments to this minimum size.
#' @param max_size Restrict fragments to this maximum size.
#' @param ... Other parameters passed to get_fragment_size
#' @return a data frame with the first three columns having the bins coordinates and
#' other columns having the WPS (raw and adjusted by coverage) and number of fragment ends (raw and adjusted by coverage).
#' @details Fragment length will extracted from the bam file according to the parameters passed to \code{\link{get_fragment_size}}, and the number of fragment ends,
#' and the Windowed Protection Score (WPS) will be computed in the binned input targets. Binning is done according to the window_size and step_size parameters.
#'
#' WPS is defined as the number of fragments completely spanning a window (bin) minus
#' the number of fragments with an endpoint within the same window as reported by Snyder et al., Cell 2016.
#'
#' The output include both the fragment end counts and the WPS in their raw format as well as after adjustment by coverage in the bin.
#'
#' Minimum and maximum bounds of the fragment size are applied before computing WPS and fragment ends counts.
#' @seealso  \code{\link{get_fragment_size}} \code{\link{bin_fragment_size}} \code{\link{summarize_fragment_size}}
#' @export
#' @importFrom rlang .data
#' @importFrom magrittr %>%
#' @examples
#' \donttest{
#' data("targets", package = "ctDNAtools")
#' bamN1 <- system.file("extdata", "N1.bam", package = "ctDNAtools")
#'
#' ## basic usage
#' analyze_fragmentation(bam = bamN1, targets = targets)
#'
#' ## more options
#' analyze_fragmentation(
#'   bam = bamN1, targets = targets,
#'   step_size = 10, window_size = 50
#' )
#' }
#'
analyze_fragmentation <- function(bam, targets, tag = "", window_size = 120,
                                  step_size = 5, min_size = 120, max_size = 180, ...) {
  assertthat::assert_that(
    assertthat::has_name(targets, c("chr", "start", "end")),
    is.data.frame(targets)
  )

  assertthat::assert_that(!missing(bam), is.character(bam), length(bam) == 1, file.exists(bam))

  assertthat::assert_that(is.character(tag), length(tag) == 1)

  assertthat::assert_that(
    is.numeric(window_size), is.numeric(step_size),
    length(window_size) == 1, length(step_size) == 1
  )

  sm <- get_bam_SM(bam = bam, tag = tag)

  ellipsis::check_dots_used()

  reads <- get_fragment_size(
    bam = bam, tag = tag, targets = targets,
    min_size = min_size, max_size = max_size, ...)

  reads_gr <- GenomicRanges::GRanges(reads$chr, IRanges::IRanges(reads$start, reads$end))

  regions <- purrr::pmap_dfr(
    list(targets$chr, targets$start, targets$end),
    function(x, y, z) {
      data.frame(
        chr = x, start = seq(y - window_size + 1, z + window_size - 1, by = step_size),
        stringsAsFactors = FALSE) %>%
        dplyr::mutate(end = .data$start + window_size)
    }
  )

  frag_ends_gr <- GenomicRanges::GRanges(
    rep(reads$chr, 2),
    IRanges::IRanges(c(reads$start, reads$end), c(reads$start, reads$end))
  )

  regions_gr <- GenomicRanges::GRanges(regions$chr, IRanges::IRanges(regions$start, regions$end))

  n_overlaps <- GenomicRanges::countOverlaps(regions_gr, reads_gr, type = "any")

  n_within_overlaps <- GenomicRanges::countOverlaps(regions_gr, reads_gr, type = "within")

  n_side_overlaps <- n_overlaps - n_within_overlaps

  n_frag_ends <- GenomicRanges::countOverlaps(regions_gr, frag_ends_gr)

  regions <- regions %>%
    dplyr::mutate(
      WPS = (n_within_overlaps - n_side_overlaps),
      WPS_adjusted = .data$WPS / ifelse(n_overlaps == 0, 1, n_overlaps),
      n_fragment_ends = n_frag_ends,
      n_fragment_ends_adjusted = n_frag_ends / ifelse(n_overlaps == 0, 1, n_overlaps),
      n_reads = n_overlaps
    )

  return(regions)
}
