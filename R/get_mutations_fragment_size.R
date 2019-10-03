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