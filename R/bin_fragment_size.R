#' Helper function to bin a variable

#' @param x the variable to be binned
#' @param from the minimum range
#' @param to the maximum range
#' @param by the break length
#' @return A numeric vector having counts within bins normalized by the sum of the variable

get_hist_bins <- function(x, from, to, by, normalized){
	
	assertthat::assert_that(is.numeric(x))

    breaks <- seq(from = from, to = to, by = by)

    if( breaks[length(breaks)] != to){

    	breaks <- c(breaks, to)
    
    }

    h <- hist(x, breaks = breaks, plot = F)
    
    if(normalized){
      
      return(h$counts/sum(h$counts))
    
    } else {
      
      return(h$counts)
    
    }
}


#' Helper function to sample proportion of a variable

#' @param x the variable to sample from.
#' @param proportion the fraction of the variable to be randomly sample.
#' @param seed the random seed.
#' @return a vector having randomly sampled proportion of the input variable.

sample_proportion <- function(x, proportion = 0.5, seed = 42) {

    sample.size <- round(length(x) * proportion)

    set.seed(seed)

    return(sample(x, size = sample.size))

}


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
#' @param augmented A logical, whether to create additional samples by sampling proportion of the reads.
#' @param n_augmentation Number of pseudo-samples to create.
#' @param augment_proportion the fraction of the reads to sample in each augmented sample.
#' @param seed The random seed to use when creating augmented samples.
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
    normalized = F, augmented = F, n_augmentation = 10, augment_proportion = 0.5, seed = 123, 
    isProperPair = NA, min_size = 1, max_size = 400, ignore_trimmed = T, 
    different_strands = T, simple_cigar = F) {

  assertthat::assert_that(is.logical(mutated_only), is.logical(augmented), is.logical(normalized),
  	    length(mutated_only) == 1, length(augmented) == 1, length(normalized) == 1)

  assertthat::assert_that(is.numeric(bin_size), bin_size %% 1 == 0,
        length(bin_size) == 1, bin_size > 0)

  assertthat::assert_that(is.numeric(n_augmentation), n_augmentation %% 1 == 0,
        length(n_augmentation) == 1)

  assertthat::assert_that(is.numeric(augment_proportion), augment_proportion < 1, 
  	    augment_proportion > 0, length(augment_proportion) == 1)

  assertthat::assert_that(is.numeric(seed), seed %% 1 == 0,
        length(seed) == 1)

  if(mutated_only) {

    assertthat::assert_that(!is.null(mutations), 
    	msg = "mutations should be specified when mutated_only is TRUE")
  
    frag_length <- get_fragment_size(bam = bam, mutations = mutations, tag = tag, 
    	isProperPair = isProperPair, min_size = min_size, max_size = max_size, 
    	ignore_trimmed = ignore_trimmed, different_strands = different_strands, 
    	simple_cigar = simple_cigar) %>%
      dplyr::filter(.data$category == "mutated")
  
  } else {

  	frag_length <- get_fragment_size(bam = bam, tag = tag, 
    	isProperPair = isProperPair, min_size = min_size, max_size = max_size, 
    	ignore_trimmed = ignore_trimmed, different_strands = different_strands, 
    	simple_cigar = simple_cigar)
  
  }
  
  message(sprintf("binning %s reads ...", nrow(frag_length)))

  out <- data.frame(counts = get_hist_bins(frag_length$size, 
     	   from = min_size, to = max_size, by = bin_size,
         normalized = normalized))

  colnames(out) <- frag_length$Sample[1]

  if(augmented) {

    set.seed(seed)

    seeds <- round(runif(n = n_augmentation, min = 0, max = 1e+08))
 
    augs <- purrr::map_dfc(seeds, 
    	~ get_hist_bins( sample_proportion(frag_length$size, proportion = augment_proportion, seed = .x),
    		from = min_size, to = max_size, by = bin_size, normalized = normalized)) %>%
       as.data.frame()

    colnames(augs) <- paste0(frag_length$Sample[1], "_resample", c(1:ncol(augs)))   

    out <- dplyr::bind_cols(out, augs) %>%
       as.data.frame()
  
  }

  return(out)

}