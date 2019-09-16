#' Counts ref and alt reads for a set of mutations
#'
#' Counts ref and alt reads for a set of mutations in a bam file
#' @param mutations A data frame with the reporter mutations. Should have the columns CHROM, POS, REF, ALT.
#' @param bam path to bam file
#' @param tag the RG tag if the bam has more than one sample
#' @param min_base_quality minimum base quality for a read to be counted
#' @param max_depth maximum depth above which sampling will happen
#' @param include_indels whether to include indels in the pileup
#' @param min_mapq the minimum mapping quality for a read to be counted
#' @return a named list contains: ref, vector of read counts of the reference alleles, and
#'         alt, vector of read counts of the alternative allele

get_mutations_read_counts <- function(mutations, bam, tag = "", min_base_quality = 20, 
            max_depth = 100000, include_indels = F, min_mapq = 30) {
  
  assertthat::assert_that(is.data.frame(mutations), assertthat::not_empty(mutations), 
                          assertthat::has_name(mutations, c("CHROM", "POS", "REF", "ALT")))
  
  ref <- purrr::pmap_dbl(list(mutations$CHROM, mutations$POS, mutations$REF),
           function(x, y, z) getReadCounts(chr = x, pos = y, base = z, bam = bam,
                                           tag = tag, min_base_quality = min_base_quality,
                                           min_mapq = min_mapq, max_depth = max_depth, include_indels = include_indels))
  
  alt <- purrr::pmap_dbl(list(mutations$CHROM, mutations$POS, mutations$ALT),
           function(x, y, z) getReadCounts(chr = x, pos = y, base = z, bam = bam,
                                           tag = tag, min_base_quality = min_base_quality,
                                           min_mapq = min_mapq, max_depth = max_depth, include_indels = include_indels))
 return(list(ref = ref, alt = alt))  
}