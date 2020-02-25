#' Computes the background mismatch rate from a bam file
#'
#' Runs through the target regions base by base counting the mismatches. Then it divides sum(mismatches)/sum(depths) for all bases in the targets
#' @param bam path to bam file
#' @param targets a data frame with the target regions. Must have three columns: chr, start and end
#' @param reference the reference genome in BSgenome format
#' @param vaf_threshold the bases with higher than this VAF threshold will be ignored in the calculation (real mutations)
#' @param tag the RG tag if the bam has more than one sample
#' @param black_list a character vector of genomic loci of format chr_pos if substitution_specific is false, or chr_pos_ref_alt if substitution_specific is true.
#' The background will be computed on the target regions after excluding black_list loci.
#' @param min_base_quality minimum base quality for a read to be counted
#' @param max_depth maximum depth above which sampling will happen
#' @param min_mapq the minimum mapping quality for a read to be counted
#' @param substitution_specific logical, whether to have the loci of black_list by substitutions.

#' @return a list containing the general mismatch rate and substitution-specific rates
#' @seealso \code{\link{create_black_list}} \code{\link{test_ctDNA}} \code{\link{create_background_panel}}
#' @details Computes the background rate of the input bam file for all bases in the specified targets. Substitutions-specific rates are also calculated.
#'
#' Genomic positions having non-reference allele frequency higher than vaf_threshold will be excluded (to exclude SNPs and real mutations).
#'
#' If a black_list is specified, the positions in the black_list (whether substitution_specific or not) will be excluded before computing the background rate.
#' @export
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @examples
#' \donttest{
#' ## Load example data
#' data("targets", package = "ctDNAtools")
#' bamT1 <- system.file("extdata", "T1.bam", package = "ctDNAtools")
#' bamN1 <- system.file("extdata", "N1.bam", package = "ctDNAtools")
#' bamN2 <- system.file("extdata", "N2.bam", package = "ctDNAtools")
#' bamN3 <- system.file("extdata", "N3.bam", package = "ctDNAtools")
#'
#' ## Use human reference genome from BSgenome.Hsapiens.UCSC.hg19 library
#' suppressMessages(library(BSgenome.Hsapiens.UCSC.hg19))
#'
#' ## basic usage
#' get_background_rate(bamT1, targets, BSgenome.Hsapiens.UCSC.hg19)
#'
#' ## more options
#' get_background_rate(bamT1, targets, BSgenome.Hsapiens.UCSC.hg19,
#'   min_base_quality = 30, min_mapq = 40, vaf_threshold = 0.05
#' )
#'
#' ## with blacklist
#' bg_panel <- create_background_panel(
#'   bam_list = c(bamN1, bamN2, bamN3),
#'   targets = targets, reference = BSgenome.Hsapiens.UCSC.hg19,
#'   substitution_specific = TRUE
#' )
#'
#' bl2 <- create_black_list(bg_panel,
#'   mean_vaf_quantile = 0.99,
#'   min_samples_one_read = 2, min_samples_two_reads = 1
#' )
#'
#' get_background_rate(bamT1, targets, BSgenome.Hsapiens.UCSC.hg19,
#'   black_list = bl2
#' )
#' }
#'
get_background_rate <- function(bam, targets, reference, vaf_threshold = 0.1, tag = "",
                                black_list = NULL, substitution_specific = TRUE, min_base_quality = 20, max_depth = 1e+05, min_mapq = 30) {
  assertthat::assert_that(class(reference) == "BSgenome")

  assertthat::assert_that(
    is.data.frame(targets), assertthat::not_empty(targets),
    assertthat::has_name(targets, c("chr", "start", "end"))
  )

  if (!is.null(black_list)) {
    assertthat::assert_that(is.character(black_list))

    if (substitution_specific) {
      assertthat::assert_that(all(purrr::map_dbl(strsplit(black_list, "_"), length) == 4),
        all(purrr::map_chr(strsplit(black_list, "_"), 3) %in% c("C", "A", "T", "G")),
        all(purrr::map_chr(strsplit(black_list, "_"), 4) %in% c("C", "A", "T", "G")),
        msg = "black_list should have characters in the format chr_pos_ref_alt when substitution_specific is true"
      )
    } else {
      assertthat::assert_that(all(purrr::map_dbl(strsplit(black_list, "_"), length) == 2),
        msg = "black_list should have characters in the format chr_pos"
      )
    }

    assertthat::assert_that(all(purrr::map_chr(strsplit(black_list, "_"), 1) %in% GenomeInfoDb::seqnames(reference)),
      msg = "Chromosomes of black_list are not in reference"
    )
  }

  gr <- GenomicRanges::reduce(GenomicRanges::GRanges(
    targets$chr,
    IRanges::IRanges(targets$start, targets$end)
  ))

  if (sum(BiocGenerics::width(gr)) < 100) {
    warning("Targets too small", immediate. = TRUE)
  }

  if (tag == "") {
    sbp <- Rsamtools::ScanBamParam(which = gr)
  } else {
    sbp <- Rsamtools::ScanBamParam(which = gr, tagFilter = list(RG = tag))
  }

  pileup_param <- Rsamtools::PileupParam(
    max_depth = max_depth, min_base_quality = min_base_quality,
    min_mapq = min_mapq, distinguish_strands = FALSE, include_deletions = FALSE, include_insertions = FALSE
  )

  p <- Rsamtools::pileup(bam, scanBamParam = sbp, pileupParam = pileup_param) %>%
    tidyr::pivot_wider(
      names_from = .data$nucleotide, values_from = .data$count,
      values_fill = list(count = 0)
    ) %>%
    as.data.frame() %>%
    dplyr::mutate(ref = as.character(BSgenome::getSeq(
      reference,
      GenomicRanges::GRanges(.data$seqnames, IRanges::IRanges(.data$pos, .data$pos))
    ))) %>%
    dplyr::mutate(depth = .data$A + .data$C + .data$G + .data$T)

  p_ann <- dplyr::mutate(p,
    refCount = purrr::map2_dbl(c(1:nrow(p)), p$ref, ~ p[.x, .y]),
    nonRefCount = .data$depth - .data$refCount
  ) %>%
    dplyr::filter((.data$nonRefCount / .data$depth) < vaf_threshold)

  if (!is.null(black_list)) {

    if (substitution_specific) {

      p <- p %>%
        tidyr::pivot_longer(names_to = "alt", cols = c("C", "A", "G", "T"), values_to = "count") %>%
        dplyr::mutate(locus = paste(.data$seqnames, .data$pos, .data$ref, .data$alt, sep = "_")) %>%
        dplyr::filter(!.data$locus %in% black_list) %>%
        dplyr::select(-.data$locus) %>%
        tidyr::pivot_wider(names_from = .data$alt, values_from = .data$count, values_fill = list(count = 0)) %>%
        dplyr::mutate(depth = .data$A + .data$C + .data$G + .data$T) %>%
        as.data.frame()

      p_ann <- dplyr::mutate(p,
        refCount = purrr::map2_dbl(c(1:nrow(p)), p$ref, ~ p[.x, .y]),
        nonRefCount = .data$depth - .data$refCount
      ) %>%
        dplyr::filter((.data$nonRefCount / .data$depth) < vaf_threshold)

    } else {
      p_ann <- p_ann %>%
        dplyr::filter(!paste(.data$seqnames, .data$pos, sep = "_") %in% black_list)
    }
  }

  rate_by_sub <- p_ann %>%
    dplyr::group_by(.data$ref) %>%
    dplyr::summarize(
      depth = sum(as.numeric(.data$depth)),
      A = sum(as.numeric(.data$A)), C = sum(as.numeric(.data$C)), G = sum(as.numeric(.data$G)),
      T = sum(as.numeric(.data$T))
    ) %>%
    as.data.frame()

  CA <- (rate_by_sub[rate_by_sub$ref == "C", "A"] + rate_by_sub[rate_by_sub$ref == "G", "T"]) /
    (rate_by_sub[rate_by_sub$ref == "C", "depth"] + rate_by_sub[rate_by_sub$ref == "G", "depth"])

  CG <- (rate_by_sub[rate_by_sub$ref == "C", "G"] + rate_by_sub[rate_by_sub$ref == "G", "C"]) /
    (rate_by_sub[rate_by_sub$ref == "C", "depth"] + rate_by_sub[rate_by_sub$ref == "G", "depth"])

  CT <- (rate_by_sub[rate_by_sub$ref == "C", "T"] + rate_by_sub[rate_by_sub$ref == "G", "A"]) /
    (rate_by_sub[rate_by_sub$ref == "C", "depth"] + rate_by_sub[rate_by_sub$ref == "G", "depth"])

  TA <- (rate_by_sub[rate_by_sub$ref == "T", "A"] + rate_by_sub[rate_by_sub$ref == "A", "T"]) /
    (rate_by_sub[rate_by_sub$ref == "T", "depth"] + rate_by_sub[rate_by_sub$ref == "A", "depth"])

  TC <- (rate_by_sub[rate_by_sub$ref == "T", "C"] + rate_by_sub[rate_by_sub$ref == "A", "G"]) /
    (rate_by_sub[rate_by_sub$ref == "T", "depth"] + rate_by_sub[rate_by_sub$ref == "A", "depth"])

  TG <- (rate_by_sub[rate_by_sub$ref == "T", "G"] + rate_by_sub[rate_by_sub$ref == "A", "C"]) /
    (rate_by_sub[rate_by_sub$ref == "T", "depth"] + rate_by_sub[rate_by_sub$ref == "A", "depth"])

  rate <- sum(as.numeric(p_ann$nonRefCount)) / sum(as.numeric(p_ann$depth))

  return(list(rate = rate, CA = CA, CG = CG, CT = CT, TA = TA, TC = TC, TG = TG))
}
