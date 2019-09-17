#' Counts ref and alt reads for a set of mutations
#'
#' Counts ref and alt reads for a set of mutations in a bam file
#' @param mutations A data frame with the reporter mutations. Should have the columns CHROM, POS, REF, ALT.
#' @param bam path to bam file
#' @param tag the RG tag if the bam has more than one sample
#' @param min_base_quality minimum base quality for a read to be counted
#' @param max_depth maximum depth above which sampling will happen
#' @param min_mapq the minimum mapping quality for a read to be counted
#' @return a named list contains: ref, vector of read counts of the reference alleles, and
#'         alt, vector of read counts of the alternative allele

get_mutations_read_counts <- function(mutations, bam, tag = "", min_base_quality = 20, 
            max_depth = 100000, min_mapq = 30) {
  
  assertthat::assert_that(is.data.frame(mutations), assertthat::not_empty(mutations), 
                          assertthat::has_name(mutations, c("CHROM", "POS", "REF", "ALT")))
  
  readcounts <- purrr::pmap_dfr(list(mutations$CHROM, mutations$POS, mutations$REF, mutations$ALT),
           function(chr, pos, ref, alt) {
              counts <- get_read_counts(chr = chr, pos = pos, bam = bam,
                         tag = tag, min_base_quality = min_base_quality,
                         min_mapq = min_mapq, max_depth = max_depth)
              data.frame(ref = counts[[ref]], alt = counts[[alt]])
            })
  
 return(list(ref = readcounts$ref, alt = readcounts$alt))  
}