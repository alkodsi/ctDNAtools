#' Filter a set of mutations 
#'
#' Filter a set of mutations given the prevalence of their mismatches in a set of bam files. Mutations that 
#' have more than min_alt_reads in more than min_samples will be removed.
#' @param mutations A data frame with the reporter mutations. Should have the columns CHROM, POS, REF, ALT.
#' @param bams a vector of paths to bam files
#' @param tags a vector of the RG tags if the bam has more than one sample
#' @param min_alt_reads the threshold of read counts showing alternative allele for a sample to be counted
#' @param min_samples the threshold of number of samples above which the mutations is filtered
#' @param min_base_quality minimum base quality for a read to be counted
#' @param max_depth maximum depth above which sampling will happen
#' @param include_indels whether to include indels in the pileup
#' @param min_mapq the minimum mapping quality for a read to be counted
#' @return a named list contains: ref, vector of read counts of the reference alleles, and
#'         alt, vector of read counts of the alternative allele
#' @export

filter_mutations <- function(mutations, bams, tags = rep("", length(bams)),
              min_alt_reads = 2, min_samples = 2, 
              min_base_quality = 20, max_depth = 100000, include_indels = F, min_mapq = 30) {
  message("Filtering mutations ...\n") 
  altMatrix <- purrr::map2_dfc(bams, tags, 
                      ~ get_mutations_read_counts(mutations = mutations, bam = .x, tag = .y,
                               min_base_quality = min_base_quality, min_mapq = min_mapq, 
                               max_depth = max_depth, include_indels = include_indels)$alt)
  
  idx <- rowSums(altMatrix > min_alt_reads) > min_samples
  message(paste("Dropped", sum(idx), "mutations\n"))
  return(mutations[!idx, ])
}