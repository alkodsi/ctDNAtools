
#' Get histogram of fragment lengths from a bam file
#'
#' The function first extracts fragment length from a bam file then computes the histogram over defined bins. If normalized is TRUE, 
#' the counts per bin will be normalized to the total read counts. Optionally,
#' it can computes the histogram of fragment lengths only for mutated reads (confirmed ctDNA molecules).

#' @param bam path to bam file.
#' @param mutations An optional data frame with mutations. Must have the columns CHROM, POS, REF, ALT.
#' @param tag the RG tag if the bam has more than one samplee.
#' @param bin_size the width of the bin (breaks) of the histogram.
#' @param custom_bins A numeric vector for custom breaks to bin the histogram of fragment length. Over-rides bin_size.
#' @param mutated_only A logical, whether to return the counts for only mutated reads. The 'mutations' input should be given when TRUE.
#' @param normalized A logical, whether to normalize the counts to the total number of reads.
#' @param min_size Integer with the lowest fragment length.
#' @param max_size Integer with the highest fragment length.
#' @param ... Other parameters passed to get_fragment_size.
#' @return A data frame with one column corresponding to the sample name from the bam. 
#' Each row has the count of fragment lengths within the bin normalized by the total number of reads.
#' @details Fragment length will extracted from the bam file according to the parameters passed to \code{\link{get_fragment_size}}, and histogram counts (optionally normalized to total counts)
#' are computed. Both equal histogram bins via bin_size and manually customized bins via custom_bins are supported. 
#'
#' By setting mutated_only to true and using an input mutations, the function would bin only fragments that support mutation variant alleles.
#' @seealso \code{\link{get_fragment_size}} \code{\link{analyze_fragmentation}} \code{\link{summarize_fragment_size}}
#' @importFrom rlang .data
#' @importFrom magrittr %>%
#' @export
#' @examples
#' \dontrun{
#' data('targets',package = 'ctDNAtools')
#' data('mutations',package = 'ctDNAtools')
#' bamT1 <- system.file('extdata', 'T1.bam', package = 'ctDNAtools')
#' 
#' ## basic usage
#' bin_fragment_size(bam = bamT1)
#' 
#' ## with normalization
#' bin_fragment_size(bam = bamT1, normalized = TRUE)
#' 
#' 
#' ## binning only mutated reads
#' bin_fragment_size(bam = bamT1, mutations = mutations, mutated_only = TRUE)
#' }

bin_fragment_size <- function(bam, mutations = NULL, tag = "", bin_size = 2, custom_bins = NULL, mutated_only = F, 
    normalized = F, min_size = 1, max_size = 400, ...) {
    
    assertthat::assert_that(is.logical(mutated_only),  
        is.logical(normalized), length(mutated_only) == 1, 
        length(normalized) == 1)
    
    assertthat::assert_that(is.numeric(bin_size), bin_size%%1 == 0, 
        length(bin_size) == 1, bin_size > 0)
    
    if(!is.null(custom_bins)){
    
        assertthat::assert_that(is.numeric(custom_bins))
    
    }

    ellipsis::check_dots_used()

    if (mutated_only) {
        
        assertthat::assert_that(!is.null(mutations), 
            msg = "mutations should be specified when mutated_only is TRUE")
        
        frag_length <- get_fragment_size(bam = bam, mutations = mutations, tag = tag, 
            min_size = min_size, max_size = max_size, ...) %>% 
        dplyr::filter(.data$category == "mutated")
        
    } else {
        
        frag_length <- get_fragment_size(bam = bam, tag = tag, min_size = min_size,
            max_size = max_size, ...)
        
    }
    
    message(sprintf("binning %s reads ...", nrow(frag_length)))
    
    histogram <- get_hist_bins(frag_length$size, from = min_size, to = max_size, 
        by = bin_size, normalized = normalized, custom_bins = custom_bins)

    out <- data.frame(Breaks = histogram$breaks, counts = histogram$counts, stringsAsFactors = F)
    
    colnames(out)[[2]] <- frag_length$Sample[[1]]
    
    return(out)
    
}
