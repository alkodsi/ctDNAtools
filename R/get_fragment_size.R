#' Get fragment lengths from bam
#'
#' A function to extract fragment lengths from a bam file. Optionally, given a mutation data frame, it can categorize read lengths 
#' in mutated vs non-mutated reads.
#' @param bam path to bam file
#' @param mutations An optional data frame with mutations. Must have the columns CHROM, POS, REF, ALT.
#' @param tag the RG tag if the bam has more than one samplee
#' @param isProperPair a logical wheter to return only proper pairs (T), only improper pairs (F), or it does not matter (NA)
#' @param low_bound Integer with the lowest fragment length
#' @param high_bound Integer with the highest fragment length
#' @return A data frame with the columns Sample (SM tag in bam), ID (read ID), size (fragment size), and category (only if mutations is provided)
#' @export 
#' @importFrom rlang .data

get_fragment_size <- function(bam, mutations = NULL, tag = "", isProperPair = NA, low_bound = 1, high_bound = 400) {

    assertthat::assert_that(!missing(bam), is.character(bam), length(bam) == 1, file.exists(bam))

    assertthat::assert_that(is.character(tag), length(tag) == 1)
    
    assertthat::assert_that(is.logical(isProperPair), is.numeric(low_bound), is.numeric(high_bound),
    	length(isProperPair) == 1, length(low_bound) == 1, length(high_bound) == 1)

    if(!is.null(mutations)) {
        
        assertthat::assert_that(is.data.frame(mutations), assertthat::not_empty(mutations), 
        
        assertthat::has_name(mutations, c("CHROM", "POS", "REF", "ALT")))

        assertthat::assert_that(all(nchar(mutations$REF) == 1),
            all(nchar(mutations$ALT) == 1),
            msg = "Only SNVs are supported") 
    
        assertthat::assert_that(is.character(mutations$REF), is.character(mutations$ALT),
            all(mutations$REF %in% c("A", "C", "T", "G")), all(mutations$ALT %in% c("A", "C", "T", "G")),
            msg = "REF and ALT in mutations should be characters having basepairs")

        assertthat::assert_that(is.numeric(mutations$POS), all(mutations$POS > 0))

        assertthat::assert_that(all(mutations$CHROM %in% get_bam_chr(bam)))
    }

    flag <- Rsamtools::scanBamFlag(isPaired = T, isProperPair = isProperPair, isUnmappedQuery = F,
    	      hasUnmappedMate = F, isFirstMateRead = T, isSecondaryAlignment = F,
    	      isNotPassingQualityControls = F, isDuplicate = F, isSupplementaryAlignment = F)
    
    if(tag != ""){
	
	    assertthat::assert_that(verify_tag(bam = bam, tag = tag), 
           msg = "Specified tag not found")
	    
	    sbp <- Rsamtools::ScanBamParam(flag = flag, 
	       tagFilter = list(RG = tag), what = Rsamtools::scanBamWhat())
    
    } else {

        sbp <- Rsamtools::ScanBamParam(flag = flag, what = Rsamtools::scanBamWhat())
    
    }
    
    sm <- get_bam_SM(bam = bam, tag = tag)
    
    scanned_bam <- Rsamtools::scanBam(bam, param = sbp)[[1]]
   
    fragment_lengths <- data.frame(Sample = sm, 
    	ID = paste(sm, scanned_bam$qname, sep = "_"), 
    	size = abs(scanned_bam$isize),
    	stringsAsFactors = F) %>%
      dplyr::filter(.data$size >= low_bound & .data$size <= high_bound)

    if(!is.null(mutations)) {
     
        read_names <- unique(unlist(get_mutations_read_names(bam = bam, tag = tag, mutations = mutations)))
        
        fragment_lengths <- dplyr::mutate(fragment_lengths, 
        	category = ifelse(.data$ID %in% read_names, "mutated", "other"))
    }

    return(fragment_lengths)
}

#' Plot density of variable
#' @param df data frame in long format
#' @param var the column of which to plot the histogram
#' @param color the column containing the categories
#' @param binwidth bin width of the histogram
#' @return a ggplot object
#' @export

plot_density <- function(df, var = "size", color = "Sample", binwidth = 2.5) {
 
    ggplot2::ggplot() + 
    ggplot2::geom_freqpoly(data = df, 
  	   ggplot2::aes_string(x = "size", y = "..density..", color = color), 
  	   binwidth = binwidth) +
    ggplot2::theme_minimal()

}



