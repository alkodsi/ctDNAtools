#' Filters a set of mutations
#'
#' This function Filters a set of mutations given the input black list or the prevalence of their mismatches in a set of bam files. Mutations that
#' have more than min_alt_reads in more than min_samples will be removed when no black list is given.
#' @param mutations A data frame with the reporter mutations. Should have the columns CHROM, POS, REF, ALT.
#' @param bams a vector of paths to bam files
#' @param black_list a character vector of genomic loci of format chr_pos to filter. If not given, the bams
#' will be scanned for mismatches in the mutations loci and the specified thresholds will be applied for filtering.
#' @param tags a vector of the RG tags if the bam has more than one sample
#' @param min_alt_reads the threshold of read counts showing alternative allele for a sample to be counted
#' @param min_samples the threshold of number of samples above which the mutations is filtered
#' @param min_base_quality minimum base quality for a read to be counted
#' @param max_depth maximum depth above which sampling will happen
#' @param min_mapq the minimum mapping quality for a read to be counted
#' @param substitution_specific logical, whether to have the loci of black_list by substitutions.

#' @return a named list contains:
#'  \itemize{
#'   \item ref: vector of read counts of the reference alleles
#'   \item alt: vector of read counts of the alternative allele
#' }
#'
#' @seealso \code{\link{create_black_list}} \code{\link{test_ctDNA}} \code{\link{create_background_panel}}
#' @details Filter a set of mutations using one of two options:
#'
#' \describe{
#' \item{1.}{By providing a black list (recommended), which includes a vector of genomic loci chr_pos when substitution_specific is false,
#'    or chr_pos_ref_alt when substitutions_specific is true. In this mode, all mutations reported in the black list are simply removed.}
#'
#' \item{2.}{By providing a set of bam files. The function will run a similar functionality to \code{\link{create_background_panel}} and filter
#'    mutations based on the min_alt_reads and min_samples criteria.}
#' }
#'
#' This function is called internally in \code{\link{test_ctDNA}} so you likely won't need to use it yourself.
#' @export
#'
#' @examples
#' data("mutations", package = "ctDNAtools")
#' filter_mutations(mutations, black_list = "chr14_106327474_C_G")
filter_mutations <- function(mutations, bams = NULL, black_list = NULL,
                             tags = rep("", length(bams)), min_alt_reads = 2,
                             min_samples = 2, min_base_quality = 20, max_depth = 1e+05, min_mapq = 30, substitution_specific = TRUE) {
  assertthat::assert_that(!is.null(bams) || !is.null(black_list))

  message("Filtering mutations ...")

  if (!is.null(black_list)) {
    assertthat::assert_that(is.character(black_list))

    if (substitution_specific) {
      assertthat::assert_that(all(purrr::map_dbl(strsplit(black_list, "_"), length) == 4),
        all(purrr::map_chr(strsplit(black_list, "_"), 3) %in% c("C", "A", "T", "G")),
        all(purrr::map_chr(strsplit(black_list, "_"), 4) %in% c("C", "A", "T", "G")),
        msg = "black_list should have characters in the format chr_pos_ref_alt"
      )

      idx <- paste(mutations$CHROM, mutations$POS, mutations$REF, mutations$ALT, sep = "_") %in% black_list
    } else {
      assertthat::assert_that(all(purrr::map_dbl(strsplit(black_list, "_"), length) == 2),
        msg = "black_list should have characters in the format chr_pos"
      )

      idx <- paste(mutations$CHROM, mutations$POS, sep = "_") %in% black_list
    }

  } else {
    assertthat::assert_that(is.character(bams))

    assertthat::assert_that(length(bams) >= min_samples)

    assertthat::assert_that(length(bams) == length(tags))


    alt_matrix <- purrr::map2_dfc(bams, tags, ~ get_mutations_read_counts(
      mutations = mutations,
      bam = .x, tag = .y, min_base_quality = min_base_quality, min_mapq = min_mapq,
      max_depth = max_depth)$alt)

    idx <- rowSums(alt_matrix > min_alt_reads) > min_samples
  }

  message(paste("Dropped", sum(idx), "mutations"))

  return(mutations[!idx, ])
}
