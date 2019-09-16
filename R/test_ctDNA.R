#' Tests the ctDNA positivity of a sample
#'
#' Given a set of reporter mutation, this functions counts the reads matching the reporter mutations in the sample to be tested, 
#' estimates the mismatch rate for the sample to be tested, and then runs a permutation test to determine whether the tested sample is positive or negative.
#' @param mutations A data frame with the reporter mutations. Should have the columns CHROM, POS, REF, ALT.
#' @param bam path to bam file
#' @param targets a data frame with the target regions. Must have three columns: chr, start and end
#' @param reference the reference genome in BSgenome format
#' @param vafThreshold the bases with higher than this VAF threshold will be ignored in the calculation (real mutations)
#' @param tag the RG tag if the bam has more than one sample
#' @param min_base_quality minimum base quality for a read to be counted
#' @param max_depth maximum depth above which sampling will happen
#' @param include_indels whether to include indels in the pileup
#' @param min_mapq the minimum mapping quality for a read to be counted
#' @param nPermutation the number of permutations
#' @param seed the random seed
#' @return a named list contains: counts, a data frame of read counts of reference and variant alleles for the reporter mutations in the tested sample,
#'         backgroundRate, a scalar estimating background rate, and pvalue, the p-value of the test
#' @export

test_ctDNA <- function(mutations, bam, targets, reference, tag = "", vafThreshold = 0.1, 
    min_base_quality = 20, max_depth = 10000, include_indels = F, min_mapq = 30, 
    nPermutation = 10000, seed = 123) {
    
    assertthat::assert_that(is.data.frame(mutations), assertthat::not_empty(mutations), 
        assertthat::has_name(mutations, c("CHROM", "POS", "REF", "ALT")))

    assertthat::assert_that(is.data.frame(targets), assertthat::not_empty(targets), 
        assertthat::has_name(targets, c("chr", "start", "end")))

    assertthat::assert_that(class(reference) == "BSgenome")
    
    message("Estimating background rate ...\n")

    bg <- getBackgroundRate(bam = bam, targets = targets, reference = reference, 
        tag = tag, vafThreshold = vafThreshold, min_base_quality = min_base_quality, 
        max_depth = max_depth, include_indels = include_indels, min_mapq = min_mapq)

    message(paste("Background rate is", bg, "\n"))
    
    message("Getting ref and alt Counts \n")

    refAltReads <- get_mutations_read_counts(mutations = mutations, bam = bam, tag = tag,
            min_base_quality = min_base_quality, max_depth = max_depth, 
            include_indels = include_indels, min_mapq = min_mapq)

    altReads <- refAltReads$alt
    
    refReads <- refAltReads$ref

    refAlt <- data.frame(Ref = refReads, Alt = altReads)
    
    message("Running permutation test \n")
    
    posTest <- positivityTest(depths = refReads + altReads, altReads = altReads, 
        rate = bg / 3, seed = seed, nPermutation = nPermutation)

    message(paste("Pvalue = ", posTest, "\n"))

    message(paste("Sample is ctDNA", ifelse(posTest < 0.05, "positive\n", "negative\n")))

    return(list(counts = refAlt, backgroundRate = bg, pvalue = posTest))
}
