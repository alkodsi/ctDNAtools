
#' Gets histogram of fragment lengths from a bam file
#'
#' The function first extracts fragment length from a bam file then computes the histogram over defined bins. If normalized is TRUE,
#' the counts per bin will be normalized to the total read counts. Optionally,
#' it can computes the histogram of fragment lengths only for mutated reads (confirmed ctDNA molecules).

#' @param bam path to bam file.
#' @param mutations An optional data frame with mutations. Must have the columns CHROM, POS, REF, ALT.
#' @param tag the RG tag if the bam has more than one sample.
#' @param targets a data frame with the target regions to restrict the reads in the bam. Must have three columns: chr, start and end
#' @param bin_size the width of the bin (breaks) of the histogram.
#' @param custom_bins A numeric vector for custom breaks to bin the histogram of fragment length. Over-rides bin_size.
#' @param normalized A logical, whether to normalize the counts to the total number of reads.
#' @param min_size Integer with the lowest fragment length.
#' @param max_size Integer with the highest fragment length.
#' @param ... Other parameters passed to get_fragment_size.
#' @return A data frame with one column for the used breaks and one having the histogram (normalized) counts. If mutations is supplied, the output will have one breaks column and three columns
#' corresponding to variant allele reads, reference allele reads, and other reads.
#' Each row has the count of fragment lengths within the bin and optionally normalized by the total number of reads.
#' @details Fragment length will extracted from the bam file according to the parameters passed to \code{\link{get_fragment_size}}, and histogram counts (optionally normalized to total counts)
#' are computed. Both equal histogram bins via bin_size and manually customized bins via custom_bins are supported.
#'
#' By using an input mutations, the function will bin separately the reads that support variant alleles, reference alleles and other reads.

#' @seealso \code{\link{get_fragment_size}} \code{\link{analyze_fragmentation}} \code{\link{summarize_fragment_size}}
#' @importFrom rlang .data
#' @importFrom magrittr %>%
#' @export
#' @examples
#' \donttest{
#' data("targets", package = "ctDNAtools")
#' data("mutations", package = "ctDNAtools")
#' bamT1 <- system.file("extdata", "T1.bam", package = "ctDNAtools")
#'
#' ## basic usage
#' bin_fragment_size(bam = bamT1)
#'
#' ## with normalization
#' bin_fragment_size(bam = bamT1, normalized = TRUE)
#'
#'
#' ## binning reads categorized based on mutations ref and alt
#' bin_fragment_size(bam = bamT1, mutations = mutations)
#'
#' ## Restrict to reads into targets
#' bin_fragment_size(bam = bamT1, targets = targets)
#' }
#'
bin_fragment_size <- function(bam, mutations = NULL, targets = NULL, tag = "", bin_size = 2, custom_bins = NULL,
                              normalized = FALSE, min_size = 1, max_size = 400, ...) {
  assertthat::assert_that(is.logical(normalized), length(normalized) == 1)

  assertthat::assert_that(
    is.numeric(bin_size), bin_size %% 1 == 0,
    length(bin_size) == 1, bin_size > 0
  )

  if (!is.null(custom_bins)) {
    assertthat::assert_that(is.numeric(custom_bins))
  }

  ellipsis::check_dots_used()

  if (!is.null(mutations)) {
    frag_length <- get_fragment_size(
      bam = bam, mutations = mutations, tag = tag,
      targets = targets, min_size = min_size, max_size = max_size, ...
    )

    histograms <- purrr::map(
      c("alt", "ref", "other"),
      ~ get_hist_bins(frag_length$size[frag_length$category == .x],
        from = min_size,
        to = max_size, by = bin_size, normalized = normalized, custom_bins = custom_bins)
    )

    out <- data.frame(
      Breaks = histograms[[1]]$breaks, purrr::map_dfc(histograms, "counts"),
      stringsAsFactors = FALSE
    )

    colnames(out)[-1] <- paste(get_bam_SM(bam = bam, tag = tag), c("alt", "ref", "other"), sep = "_")
  } else {
    frag_length <- get_fragment_size(
      bam = bam, tag = tag, min_size = min_size,
      max_size = max_size, targets = targets, ...
    )

    histogram <- get_hist_bins(frag_length$size,
      from = min_size, to = max_size,
      by = bin_size, normalized = normalized, custom_bins = custom_bins
    )

    out <- data.frame(Breaks = histogram$breaks, counts = histogram$counts, stringsAsFactors = FALSE)

    colnames(out)[[2]] <- get_bam_SM(bam = bam, tag = tag)
  }

  return(out)
}
