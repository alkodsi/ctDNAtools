#' Filter a set of mutations 
#'
#' Filter a set of mutations given the input black list or the prevalence of their mismatches in a set of bam files. Mutations that 
#' have more than min_alt_reads in more than min_samples will be removed when no black list is given.
#' @param mutations A data frame with the reporter mutations. Should have the columns CHROM, POS, REF, ALT.
#' @param bams a vector of paths to bam files
#' @param black_list a character vector of genomic loci of format chr_pos to filter. If not given, the bams
#' will be scanned for mismatches in the mutations loci and the specified thresholds will be applied for filtering.
#' @param tags a vector of the RG tags if the bam has more than one sample
#' @param min_alt_reads the threshold of read counts showing alternative allele for a sample to be counted
#' @param min_samples the threshold of number of samples above which the mutations is filtered
#' @param min_base_quality minimum base quality for a read to be counted
#' @param max_depth maximum depth above which sampling will happen
#' @param min_mapq the minimum mapping quality for a read to be counted
#' @return a named list contains: ref, vector of read counts of the reference alleles, and
#'         alt, vector of read counts of the alternative allele
#' @export

filter_mutations <- function(mutations, bams = NULL, black_list = NULL, 
    tags = rep("", length(bams)), min_alt_reads = 2, 
    min_samples = 2, min_base_quality = 20, max_depth = 1e+05, min_mapq = 30) {
    
    assertthat::assert_that(!is.null(bams) || !is.null(black_list))

    message("Filtering mutations ...")

    if(!is.null(black_list)){
        
        assertthat::assert_that(is.character(black_list))
        
        assertthat::assert_that(all(purrr::map_dbl(strsplit(black_list, "_"),length) == 2),
            msg = "black_list should have characters in the format chr_pos")

        idx <- paste(mutations$CHROM, mutations$POS, sep = "_") %in% black_list

    } else {

        assertthat::assert_that(is.character(bams))

        assertthat::assert_that(length(bams) >= min_samples)
        
        assertthat::assert_that(length(bams) == length(tags))
    
    
        altMatrix <- purrr::map2_dfc(bams, tags, ~ get_mutations_read_counts(mutations = mutations, 
             bam = .x, tag = .y, min_base_quality = min_base_quality, min_mapq = min_mapq, 
             max_depth = max_depth)$alt)
    
        idx <- rowSums(altMatrix > min_alt_reads) > min_samples
    }

    message(paste("Dropped", sum(idx), "mutations"))
    
    return(mutations[!idx, ])
}
