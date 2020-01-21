
#' Creates a background panel instance from a bam file (e.g. healthy samples)
#'
#' This function scans the targets regions in one bam file, and reports the number of reference, non-reference reads for each loci
#' in addition to the non-reference (VAF) allele frequency. Loci with VAF higher than vaf_threshold are masked with NA. This function is
#' used internally by create_background_panel

#' @param bam A character specifying the path to bam file.
#' @param targets The targets data frame must have the columns chr, start and end.
#' @param reference The reference genome as BSgenome object.
#' @param vaf_threshold Loci with the fraction of non-reference reads above this value are masked with NA.
#' @param tag The RG tag in the bam file. Empty string defaults to using the whole bam.
#' @param min_base_quality The minimum base quality to count a read for a loci.
#' @param max_depth Maximum depth for the pileup
#' @param min_mapq The minimum mapping quality to count a read for a loci
#' @param substitution_specific logical, whether to have the loci by substitutions.
#' @return A named list having depth, alt and vaf data frames. Each has the same order of loci in rows and the input sample in columns.
#' @keywords internal

create_background_panel_instance <- function(bam, targets, reference, vaf_threshold = 0.05, tag = "",
                                             min_base_quality = 20, max_depth = 1e+05, min_mapq = 30, substitution_specific = TRUE) {
  gr <- GenomicRanges::reduce(GenomicRanges::GRanges(targets$chr, IRanges::IRanges(targets$start, targets$end)))

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
      GenomicRanges::GRanges(.data$seqnames, IRanges::IRanges(.data$pos, .data$pos))))
    ) %>%
    dplyr::mutate(depth = .data$A + .data$C + .data$G + .data$T)

  if (substitution_specific) {

    p_ann <- p %>%
      dplyr::mutate(
        VAFC = .data$C / .data$depth,
        VAFT = .data$T / .data$depth,
        VAFG = .data$G / .data$depth,
        VAFA = .data$A / .data$depth
      ) %>%
      tidyr::pivot_longer(names_to = "alt", cols = dplyr::starts_with("VAF"), values_to = "vaf") %>%
      dplyr::mutate(alt = gsub("VAF", "", .data$alt)) %>%
      dplyr::filter(.data$ref != .data$alt) %>%
      dplyr::mutate(
        Locus = paste(.data$seqnames, .data$pos, .data$ref, .data$alt, sep = "_"),
        nonRefCount = .data$depth * .data$vaf
      ) %>%
      dplyr::select(.data$Locus, .data$depth, .data$nonRefCount, .data$vaf)

  } else {

    p_ann <- dplyr::mutate(p,
      refCount = purrr::map2_dbl(c(1:nrow(p)), p$ref, ~ p[.x, .y]),
      nonRefCount = .data$depth - .data$refCount, vaf = .data$nonRefCount / .data$depth
    ) %>%
      dplyr::mutate(
        nonRefCount = ifelse(.data$vaf >= vaf_threshold, NA, .data$nonRefCount),
        depth = ifelse(.data$vaf >= vaf_threshold, NA, .data$depth),
        vaf = ifelse(.data$vaf >= vaf_threshold, NA, .data$vaf),
        Locus = paste(.data$seqnames, .data$pos, sep = "_")
      ) %>%
      dplyr::select(.data$Locus, .data$depth, .data$nonRefCount, .data$vaf)
  }

  out <- list(
    depth = data.frame(Locus = p_ann$Locus, depth = p_ann$depth, stringsAsFactors = FALSE),
    alt = data.frame(Locus = p_ann$Locus, alt = p_ann$nonRefCount, stringsAsFactors = FALSE),
    vaf = data.frame(Locus = p_ann$Locus, vaf = p_ann$vaf, stringsAsFactors = FALSE)
  )

  return(out)
}
