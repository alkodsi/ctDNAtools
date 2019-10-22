#' Tests the ctDNA positivity of a sample
#'
#' Given a set of reporter mutation, this functions counts the reads matching the reporter mutations in the sample to be tested, 
#' estimates the mismatch rate for the sample to be tested, and then runs a Monte Carlo simulation test to determine whether the tested sample is positive or negative.
#' @param mutations A data frame with the reporter mutations. Should have the columns CHROM, POS, REF, ALT.
#' @param bam path to bam file
#' @param targets a data frame with the target regions. Must have three columns: chr, start and end
#' @param reference the reference genome in BSgenome format
#' @param vaf_threshold the bases with higher than this VAF threshold will be ignored in the calculation (real mutations)
#' @param tag the RG tag if the bam has more than one sample
#' @param ID_column The name of the column that contains the ID of mutations in phase. All mutations in Phase should have the same ID in that column. 
#' Will lead to considerable slow down when provided.
#' @param use_unique_molecules A logical. If true, reads mapping to multiple mutations will 
#' be counted only once for the mutation that appears first in mutations data frame.
#' @param min_base_quality minimum base quality for a read to be counted
#' @param max_depth maximum depth above which sampling will happen
#' @param min_mapq the minimum mapping quality for a read to be counted
#' @param bam_list A vector containing the paths to bam files used to filter mutations. Mutations that have more than min_alt_reads in more than min_samples will be filtered.
#' @param bam_list_tags the RG tags for the bams included in bams list.
#' @param min_alt_reads When bam_list is provided, this sets the minimum number of alternative allele reads for a sample to be counted.
#' @param min_samples When number of samples having more than min_alt_reads exceeds this number, the mutation will be filtered.
#' @param by_substitution boolean whether to run the test according to substitution-specific background rate
#' @param n_simulations the number of simulations
#' @param seed the random seed
#' @return a named list contains: counts, a data frame of read counts of reference and variant alleles for the reporter mutations in the tested sample,
#'         backgroundRate, a list of substituion-specific background rate, and pvalue, the p-value of the test
#' @export

test_ctDNA <- function(mutations, bam, targets, reference, tag = "", ID_column = NULL, use_unique_molecules = T,
    vaf_threshold = 0.1, min_base_quality = 20, max_depth = 1e+05, min_mapq = 30, bam_list = character(), 
    bam_list_tags = rep("", length(bam_list)), min_alt_reads = 1, min_samples = 1, 
    by_substitution = F, n_simulations = 10000, seed = 123) {
    
    assertthat::assert_that(!missing(mutations), !missing(bam), !missing(targets), 
        !missing(reference), msg = "mutations, bam, targets and reference are all required")
    
    assertthat::assert_that(is.data.frame(mutations), assertthat::not_empty(mutations), 
        assertthat::has_name(mutations, c("CHROM", "POS", "REF", "ALT")))
    
    assertthat::assert_that(is.data.frame(targets), assertthat::not_empty(targets), 
        assertthat::has_name(targets, c("chr", "start", "end")))
    
    assertthat::assert_that(class(reference) == "BSgenome")
    
    assertthat::assert_that(all(nchar(mutations$REF) == 1), all(nchar(mutations$ALT) == 1),
        msg = "Only SNVs are supported")
    
    assertthat::assert_that(is.character(mutations$REF), is.character(mutations$ALT), 
        all(mutations$REF %in% c("A", "C", "T", "G")), all(mutations$ALT %in% c("A", "C", "T", "G")),
        msg = "REF and ALT in mutations should be characters having basepairs")
    
    assertthat::assert_that(is.numeric(mutations$POS), all(mutations$POS > 0))
    
    assertthat::assert_that(is.numeric(targets$start), all(targets$start > 0), 
        is.numeric(targets$end), all(targets$end > 0))
    
    assertthat::assert_that(all(mutations$CHROM %in% GenomeInfoDb::seqnames(reference)), 
        all(targets$chr %in% GenomeInfoDb::seqnames(reference)), 
        msg = "Chromosomes in mutations and/or targets don't match the specified reference")
    
    assertthat::assert_that(is.character(bam), length(bam) == 1, file.exists(bam))
    
    assertthat::assert_that(is.character(tag), length(tag) == 1)
    
    assertthat::assert_that(all(get_bam_chr(bam) %in% GenomeInfoDb::seqnames(reference)), 
        msg = "Chromosomes in bam file don't match the specified reference")
    
    assertthat::assert_that(is.numeric(vaf_threshold), vaf_threshold > 0, 
        vaf_threshold <= 1, length(vaf_threshold) == 1)
    
    assertthat::assert_that(is.numeric(min_base_quality), length(min_base_quality) == 1,
        min_base_quality%%1 == 0, min_base_quality > 0)
    
    assertthat::assert_that(is.numeric(max_depth), length(max_depth) == 1, 
        max_depth%%1 == 0, max_depth > 0)
    
    assertthat::assert_that(is.numeric(min_mapq), length(min_mapq) == 1, 
        min_mapq%%1 == 0, min_mapq > 0)
    
    assertthat::assert_that(is.numeric(n_simulations), length(n_simulations) == 1, 
        n_simulations%%1 == 0, n_simulations > 0)
    
    assertthat::assert_that(is.numeric(seed), length(seed) == 1, seed%%1 == 0, 
        seed >= 0)
    
    assertthat::assert_that(is.logical(by_substitution), length(by_substitution) == 1)
    
    if (tag != "") {
        
        assertthat::assert_that(verify_tag(bam = bam, tag = tag), msg = "Specified tag not found")
    
    }
    
     if (!is.null(ID_column)) {
        
        assertthat::assert_that(assertthat::has_name(mutations, ID_column))
        
        assertthat::assert_that(is.character(mutations[,ID_column]))
    
    }

    if (length(bam_list) > 0) {
        
        assertthat::assert_that(is.numeric(min_alt_reads), length(min_alt_reads) == 1,
            min_alt_reads%%1 == 0, min_alt_reads >= 0)
        
        assertthat::assert_that(is.numeric(min_samples), length(min_samples) == 1, 
            min_samples%%1 == 0, min_samples >= 0)
        
        assertthat::assert_that(is.character(bam_list), all(file.exists(bam_list)))
        
        bam_list_chr <- purrr::map(bam_list, get_bam_chr)
        
        assertthat::assert_that(all(purrr::map_lgl(bam_list_chr, ~ all(.x %in% GenomeInfoDb::seqnames(reference)))), 
            msg = "The chromosomes in at least one of the specified bams in bam_list don't match the reference")
        
        assertthat::assert_that(is.character(bam_list_tags), 
            length(bam_list_tags) == length(bam_list))
        
        if (any(bam_list_tags != "")) {
            
            tag_verification <- purrr::map2_lgl(bam_list[bam_list_tags != ""],
                 bam_list_tags[bam_list_tags != ""], verify_tag)
            
            assertthat::assert_that(all(tag_verification), 
                msg = paste("specified tag for", paste(bam_list[bam_list_tags != ""][!tag_verification], collapse = " , "), 
                "is not correct"))
        
        }
    }
    
    sm <- get_bam_SM(bam = bam, tag = tag)
    
    message(paste("Analyzing Sample", sm, "..."))
    
    if (length(bam_list) > 0) {
        
        mutations <- filter_mutations(mutations = mutations, bams = bam_list, tags = bam_list_tags, 
            min_alt_reads = min_alt_reads, min_samples = min_samples, min_base_quality = min_base_quality, 
            max_depth = max_depth, min_mapq = min_mapq)
        
    }
    
    
    message("Estimating background rate ...")
    
    bg <- get_background_rate(bam = bam, targets = targets, reference = reference, 
        tag = tag, vaf_threshold = vaf_threshold, min_base_quality = min_base_quality, 
        max_depth = max_depth, min_mapq = min_mapq)
    
    
    
    message("Getting ref and alt Counts ...")
    
    if(!is.null(ID_column)){
        
        message("merging mutations in phase ...")

        refAltReads <- merge_mutations_in_phase(mutations = mutations, bam = bam,
            tag = tag, min_base_quality = min_base_quality, ID_column = ID_column, 
            use_unique_molecules = use_unique_molecules)

        prob_purification <- refAltReads$purification_prob
        
        refAltReads <- refAltReads$out
         
        bg$rate <- bg$rate * (1 - prob_purification)

    } else {
       
        refAltReads <- get_mutations_read_counts(mutations = mutations, bam = bam, tag = tag, 
             min_base_quality = min_base_quality, max_depth = max_depth, min_mapq = min_mapq)

    }

    
    altReads <- refAltReads$alt
    
    refReads <- refAltReads$ref
    
    refAlt <- data.frame(Ref = refReads, Alt = altReads)
    
    
    message("Running Monte Carlo simulations")
    
    if (by_substitution && is.null(ID_column)) {
        
        substitutions <- paste0(mutations$REF, mutations$ALT)
        
        posTest <- positivity_test(depths = refReads + altReads, altReads = altReads, 
            substitutions = substitutions, rate = bg, seed = seed, n_simulations = n_simulations)
        
    } else {
        
        posTest <- positivity_test(depths = refReads + altReads, altReads = altReads, 
            rate = bg, seed = seed, n_simulations = n_simulations)
        
    }
    
    message(paste("Pvalue = ", posTest))
    
    message(paste("Sample", sm, "is ctDNA", ifelse(posTest < 0.05, "positive", "negative")))
    
    return(list(counts = refAlt, backgroundRate = bg, pvalue = posTest))
}
