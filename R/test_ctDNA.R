#' Tests the ctDNA positivity of a sample
#'
#' Given a set of reporter mutation, this functions counts the reads matching the reporter mutations in the sample to be tested, 
#' estimates the mismatch rate for the sample to be tested, and then runs a permutation test to determine whether the tested sample is positive or negative.
#' @param mutations A data frame with the reporter mutations. Should have the columns CHROM, POS, REF, ALT.
#' @param bam path to bam file
#' @param targets a data frame with the target regions. Must have three columns: chr, start and end
#' @param reference the reference genome in BSgenome format
#' @param vaf_threshold the bases with higher than this VAF threshold will be ignored in the calculation (real mutations)
#' @param tag the RG tag if the bam has more than one sample
#' @param min_base_quality minimum base quality for a read to be counted
#' @param max_depth maximum depth above which sampling will happen
#' @param min_mapq the minimum mapping quality for a read to be counted
#' @param bam_list A vector containing the paths to bam files used to filter mutations. Mutations that have more than min_alt_reads in more than min_samples will be filtered.
#' @param bam_list_tags the RG tags for the bams included in bams list.
#' @param min_alt_reads When bam_list is provided, this sets the minimum number of alternative allele reads for a sample to be counted.
#' @param min_samples When number of samples having more than min_alt_reads exceeds this number, the mutation will be filtered.
#' @param by_substitution boolean whether to run the test according to substitution-specific background rate
#' @param n_permutation the number of permutations
#' @param seed the random seed
#' @return a named list contains: counts, a data frame of read counts of reference and variant alleles for the reporter mutations in the tested sample,
#'         backgroundRate, a list of substituion-specific background rate, and pvalue, the p-value of the test
#' @export

test_ctDNA <- function(mutations, bam, targets, reference, tag = "", vaf_threshold = 0.1, 
    min_base_quality = 20, max_depth = 1e+05, min_mapq = 30, 
    bam_list = character(), bam_list_tags = rep("",length(bam_list)), min_alt_reads = 1, min_samples = 1,
    by_substitution = T, n_permutation = 1e+04, seed = 123) {
    
    assertthat::assert_that(is.data.frame(mutations), assertthat::not_empty(mutations), 
        assertthat::has_name(mutations, c("CHROM", "POS", "REF", "ALT")))

    assertthat::assert_that(is.data.frame(targets), assertthat::not_empty(targets), 
        assertthat::has_name(targets, c("chr", "start", "end")))

    assertthat::assert_that(class(reference) == "BSgenome")
     
    if(length(bam_list) > 0){
       
        mutations <- filter_mutations(mutations = mutations, bams = bam_list, tags = bam_list_tags,
          min_alt_reads = min_alt_reads, min_samples = min_samples, min_base_quality = min_base_quality,
          max_depth = max_depth, min_mapq = min_mapq)
    }

    subs <- paste0(mutations$REF, mutations$ALT)
    assertthat::assert_that(all(nchar(subs) == 2), msg = "Only SNVs are supported")

    message("Estimating background rate ...\n")

    bg <- get_background_rate(bam = bam, targets = targets, reference = reference, 
        tag = tag, vaf_threshold = vaf_threshold, min_base_quality = min_base_quality, 
        max_depth = max_depth, min_mapq = min_mapq)

    
    message("Getting ref and alt Counts \n")

    refAltReads <- get_mutations_read_counts(mutations = mutations, bam = bam, tag = tag,
            min_base_quality = min_base_quality, max_depth = max_depth, min_mapq = min_mapq)

    altReads <- refAltReads$alt
    
    refReads <- refAltReads$ref

    refAlt <- data.frame(Ref = refReads, Alt = altReads)
    
    message("Running permutation test \n")
    
    if(by_substitution){
           
      substitutions <- dplyr::case_when(subs %in% c("CT", "GA") ~ "CT",
                  subs %in% c("CA", "GT") ~ "CA",
                  subs %in% c("CG", "GC") ~ "CG",
                  subs %in% c("TA", "AT") ~ "TA",
                  subs %in% c("TC", "AG") ~ "TC",
                  subs %in% c("TG", "AC") ~ "TG")

      posTest <- positivity_test(depths = refReads + altReads, altReads = altReads, 
        substitutions = substitutions, rate = bg, seed = seed, n_permutation = n_permutation)
    
    } else {
    
      posTest <- positivity_test(depths = refReads + altReads, altReads = altReads, 
        rate = bg, seed = seed, n_permutation = n_permutation)
    
    }
    
    message(paste("Pvalue = ", posTest, "\n"))

    message(paste("Sample is ctDNA", ifelse(posTest < 0.05, "positive\n", "negative\n")))

    return(list(counts = refAlt, backgroundRate = bg, pvalue = posTest))
}
