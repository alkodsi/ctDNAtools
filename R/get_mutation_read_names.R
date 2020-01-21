
#' Gets names of the reads showing reference and alternative alleles of a mutation
#'
#' Extract the names of the reads in a bam file that support the variant allele of a single mutation
#' @param bam path to bam file.
#' @param chr Chromosome name for the mutation.
#' @param pos Chromosome position for the mutation.
#' @param ref The reference allele of the mutation.
#' @param alt The alternative allele of the mutation.
#' @param tag the RG tag if the bam has more than one sample.
#' @param min_base_quality integer specifying the minimum base quality for reads to be included.
#' @param min_mapq integer specifying the minimum mapping quality for reads to be included.
#' @return A character vector having the read names
#' @importFrom rlang .data
#' @keywords internal

get_mutation_read_names <- function(bam, chr, pos, ref, alt, tag = "", min_base_quality = 20, min_mapq = 30) {
  assertthat::assert_that(!missing(bam), is.character(bam), length(bam) == 1, file.exists(bam))

  assertthat::assert_that(is.character(tag), length(tag) == 1)

  assertthat::assert_that(is.numeric(pos), pos > 0, pos %% 1 == 0, length(pos) == 1)

  assertthat::assert_that(is.character(chr), length(chr) == 1, chr %in% get_bam_chr(bam))

  assertthat::assert_that(is.character(ref), length(ref) == 1, ref %in% c("C", "T", "G", "A"))

  assertthat::assert_that(is.character(alt), length(alt) == 1, alt %in% c("C", "T", "G", "A"))

  gr <- GenomicRanges::GRanges(chr, IRanges::IRanges(pos, pos))

  if (tag == "") {
    sbp <- Rsamtools::ScanBamParam(
      flag = Rsamtools::scanBamFlag(
        isUnmappedQuery = FALSE, isSecondaryAlignment = FALSE,
        isNotPassingQualityControls = FALSE, isDuplicate = FALSE
      ),
      which = gr, mapqFilter = min_mapq
    )

    stacked_strings <- GenomicAlignments::stackStringsFromBam(bam,
      use.names = TRUE, param = sbp)

    stacked_qual <- GenomicAlignments::stackStringsFromBam(bam,
      use.names = TRUE, what = "qual", param = sbp)

  } else {
    assertthat::assert_that(verify_tag(bam = bam, tag = tag), msg = "Specified tag not found")

    sbp <- Rsamtools::ScanBamParam(
      flag = Rsamtools::scanBamFlag(
        isUnmappedQuery = FALSE, isSecondaryAlignment = FALSE,
        isNotPassingQualityControls = FALSE, isDuplicate = FALSE
      ),
      which = gr, tagFilter = list(RG = tag), mapqFilter = min_mapq
    )

    stacked_strings <- GenomicAlignments::stackStringsFromBam(bam,
      use.names = TRUE, param = sbp)

    stacked_qual <- GenomicAlignments::stackStringsFromBam(bam,
      use.names = TRUE, what = "qual", param = sbp)
  }

  sm <- get_bam_SM(bam = bam, tag = tag)

  if (length(stacked_strings) != 0) {
    out <- data.frame(
      ID = paste(sm, names(stacked_strings), sep = "_"),
      seq = as.data.frame(stacked_strings)[, 1],
      qual = as.data.frame(methods::as(Biostrings::PhredQuality(stacked_qual), "IntegerList"))$value,
      stringsAsFactors = FALSE
    ) %>%
      dplyr::filter(.data$qual >= min_base_quality)

    alt <- unique(out[out$seq == alt, "ID"])
    ref <- unique(out[out$seq == ref, "ID"])
    return(list(ref = ref, alt = alt))
  } else {
    return(list(ref = character(0), alt = character(0)))
  }
}
