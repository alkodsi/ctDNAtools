#' Gets read counts for a specific locus in the genome
#'
#' Uses samtools pileup to get the read counts for each base in the genomic position specified
#' @param chr chromosome name
#' @param pos genomic coordinate
#' @param bam path to bam file
#' @param tag the RG tag if the bam has more than one sample
#' @param min_base_quality minimum base quality for a read to be counted
#' @param max_depth maximum depth above which sampling will happen
#' @param min_mapq the minimum mapping quality for a read to be counted
#' @return a list, number of reads for each of the four basepairs
#' @keywords internal

get_read_counts <- function(chr, pos, bam, tag = "", min_base_quality = 20, max_depth = 1e+05,
                            min_mapq = 30) {
  gr <- GenomicRanges::GRanges(chr, IRanges::IRanges(pos, pos))

  if (tag == "") {
    sbp <- Rsamtools::ScanBamParam(which = gr)
  } else {
    sbp <- Rsamtools::ScanBamParam(which = gr, tagFilter = list(RG = tag))
  }

  pileup_param <- Rsamtools::PileupParam(
    max_depth = max_depth, min_base_quality = min_base_quality,
    min_mapq = min_mapq, distinguish_strands = FALSE, include_deletions = FALSE, include_insertions = FALSE
  )

  p <- Rsamtools::pileup(bam, scanBamParam = sbp, pileupParam = pileup_param)

  cbase <- ifelse("C" %in% p$nucleotide, p[p$nucleotide == "C", "count"], 0)
  gbase <- ifelse("G" %in% p$nucleotide, p[p$nucleotide == "G", "count"], 0)
  abase <- ifelse("A" %in% p$nucleotide, p[p$nucleotide == "A", "count"], 0)
  tbase <- ifelse("T" %in% p$nucleotide, p[p$nucleotide == "T", "count"], 0)

  return(list(A = abase, C = cbase, G = gbase, T = tbase))
}
