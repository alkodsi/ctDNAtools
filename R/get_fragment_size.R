#' Get fragment lengths from bam
#'
#' A function to extract fragment lengths from a bam file. Optionally, given a mutation data frame, it can categorize read lengths 
#' in mutated vs non-mutated reads.
#' @param bam path to bam file.
#' @param mutations An optional data frame with mutations. Must have the columns CHROM, POS, REF, ALT.
#' @param tag the RG tag if the bam has more than one samplee.
#' @param isProperPair a logical wheter to return only proper pairs (T), only improper pairs (F), or it does not matter (NA).
#' @param min_size Integer with the lowest fragment length.
#' @param max_size Integer with the highest fragment length.
#' @param ignore_trimmed logical, whether to remove reads that have been hard trimmed.
#' @param different_strands logical, whether to keep only reads whose mates map to different strand.
#' @param simple_cigar logical, whether to include only reads with simple cigar.
#' @return A data frame with the columns Sample (SM tag in bam), ID (read ID), size (fragment size), and category (only if mutations is provided).
#' @export 
#' @importFrom rlang .data
#' @importFrom magrittr %>%

get_fragment_size <- function(bam, mutations = NULL, tag = "", isProperPair = NA, 
	min_size = 1, max_size = 400, ignore_trimmed = T, different_strands = T, simple_cigar = F) {

    assertthat::assert_that(!missing(bam), is.character(bam), length(bam) == 1, file.exists(bam))

    assertthat::assert_that(is.character(tag), length(tag) == 1)
    
    assertthat::assert_that(is.logical(isProperPair), is.logical(ignore_trimmed), 
    	is.logical(different_strands), is.logical(simple_cigar),
    	is.numeric(min_size), is.numeric(max_size),
    	length(isProperPair) == 1, length(ignore_trimmed) == 1, length(different_strands) == 1, 
    	length(simple_cigar) == 1, length(min_size) == 1, length(max_size) == 1)

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
	    
	    sbp <- Rsamtools::ScanBamParam(flag = flag, mapqFilter = 30, simpleCigar = simple_cigar,
	       tagFilter = list(RG = tag), what = c("qname","flag","qwidth","isize"))
    
    } else {

        sbp <- Rsamtools::ScanBamParam(flag = flag, mapqFilter = 30, simpleCigar = simple_cigar,
           what = c("qname","flag","qwidth","isize"))
    
    }
    
    sm <- get_bam_SM(bam = bam, tag = tag)
    
    scanned_bam <- Rsamtools::scanBam(bam, param = sbp)[[1]] %>%
        dplyr::bind_cols() %>% as.data.frame()

    if(different_strands){
       
        flag_matrix <- Rsamtools::bamFlagAsBitMatrix(as.integer(scanned_bam$flag))
    
        scanned_bam <- scanned_bam[xor(flag_matrix[,'isMinusStrand'], flag_matrix[,'isMateMinusStrand']),]
    }

    if(ignore_trimmed){
    	
    	read_length <- max(scanned_bam$qwidth)
        
        scanned_bam <- scanned_bam %>%
          filter(.data$qwidth == read_length)
    }

    fragment_lengths <- data.frame(Sample = sm, 
    	ID = paste(sm, scanned_bam$qname, sep = "_"),
    	size = abs(scanned_bam$isize),
    	stringsAsFactors = F) %>%
      dplyr::filter(.data$size >= min_size & .data$size <= max_size)

    if(!is.null(mutations)) {
     
        read_names <- unique(unlist(get_mutations_read_names(bam = bam, tag = tag, mutations = mutations)))
        
        fragment_lengths <- dplyr::mutate(fragment_lengths, 
        	category = ifelse(.data$ID %in% read_names, "mutated", "other"))
    }

    return(fragment_lengths)
}


#' Get reads fragment lengths for a list of mutations
#'
#' The function extracts the fragment lengths for the reads holding alternative allele for each mutation
#' in the mutations data frame. 
#' @param bam path to bam file.
#' @param mutations Data frame with mutations. Must have the columns CHROM, POS, REF, ALT.
#' @param tag the RG tag if the bam has more than one samplee.
#' @param isProperPair a logical wheter to return only proper pairs (T), only improper pairs (F), or it does not matter (NA).
#' @param min_size Integer with the lowest fragment length.
#' @param max_size Integer with the highest fragment length.
#' @param ignore_trimmed logical, whether to remove reads that have been hard trimmed.
#' @param different_strands logical, whether to keep only reads whose mates map to different strand.
#' @param simple_cigar logical, whether to include only reads with simple cigar.
#' @return A list with length equal to the number of mutations. Each element contains an integer vector of fragment lengths
#' @export 

get_mutations_fragment_size <- function(bam, mutations, tag = "", isProperPair = NA, 
	min_size = 1, max_size = 400, ignore_trimmed = T, different_strands = T, simple_cigar = F) {

  assertthat::assert_that(!missing(mutations), !missing(bam))

  frag_size <- get_fragment_size(bam = bam, tag = tag, isProperPair = isProperPair,
  	    min_size = min_size, max_size = max_size, ignore_trimmed = ignore_trimmed, 
  	    different_strands = different_strands, simple_cigar = simple_cigar)

  read_names <- get_mutations_read_names(bam = bam, mutations = mutations, tag = tag)

  mutation_frag_size <- purrr::map(read_names, ~ frag_size[ frag_size$ID %in% .x, "size"])

  return(mutation_frag_size)
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
  	   ggplot2::aes_string(x = var, y = "..density..", color = color), 
  	   binwidth = binwidth) +
    ggplot2::theme_minimal()

}