
#' Get histogram of fragment lengths from a bam file
#'
#' The function first extracts fragment length from a bam file then computes the histogram over defined bins. If normalized is TRUE, 
#' the counts per bin will be normalized to the total read counts. Optionally,
#' it can computes the histogram of fragment lengths only for mutated reads (confirmed ctDNA molecules).

#' @param bam path to bam file.
#' @param mutations An optional data frame with mutations. Must have the columns CHROM, POS, REF, ALT.
#' @param tag the RG tag if the bam has more than one samplee.
#' @param bin_size the width of the bin (breaks) of the histogram.
#' @param mutated_only A logical, whether to return the counts for only mutated reads. The 'mutations' input should be given when TRUE.
#' @param normalized A logical, whether to normalize the counts to the total number of reads.
#' @param isProperPair a logical wheter to return only proper pairs (T), only improper pairs (F), or it does not matter (NA).
#' @param min_size Integer with the lowest fragment length.
#' @param max_size Integer with the highest fragment length.
#' @param ignore_trimmed logical, whether to remove reads that have been hard trimmed.
#' @param different_strands logical, whether to keep only reads whose mates map to different strand.
#' @param simple_cigar logical, whether to include only reads with simple cigar.
#' @return A data frame with one column corresponding to the sample name from the bam. 
#' Each row has the count of fragment lengths within the bin normalized by the total number of reads.
#' @importFrom rlang .data
#' @importFrom magrittr %>%
#' @export

bin_fragment_size <- function(bam, mutations = NULL, tag = "", bin_size = 2, mutated_only = F, 
    normalized = F, isProperPair = NA, min_size = 1, max_size = 400, ignore_trimmed = T, 
    different_strands = T, simple_cigar = F) {
    
    assertthat::assert_that(is.logical(mutated_only),  
        is.logical(normalized), length(mutated_only) == 1, 
        length(normalized) == 1)
    
    assertthat::assert_that(is.numeric(bin_size), bin_size%%1 == 0, 
        length(bin_size) == 1, bin_size > 0)
    
    if (mutated_only) {
        
        assertthat::assert_that(!is.null(mutations), 
            msg = "mutations should be specified when mutated_only is TRUE")
        
        frag_length <- get_fragment_size(bam = bam, mutations = mutations, tag = tag, 
            isProperPair = isProperPair, min_size = min_size, max_size = max_size, 
            ignore_trimmed = ignore_trimmed, different_strands = different_strands, 
            simple_cigar = simple_cigar) %>% dplyr::filter(.data$category == "mutated")
        
    } else {
        
        frag_length <- get_fragment_size(bam = bam, tag = tag, isProperPair = isProperPair, 
            min_size = min_size, max_size = max_size, ignore_trimmed = ignore_trimmed, 
            different_strands = different_strands, simple_cigar = simple_cigar)
        
    }
    
    message(sprintf("binning %s reads ...", nrow(frag_length)))
    
    out <- data.frame(counts = get_hist_bins(frag_length$size, from = min_size, to = max_size, 
        by = bin_size, normalized = normalized))
    
    colnames(out) <- frag_length$Sample[1]
    
    return(out)
    
}
