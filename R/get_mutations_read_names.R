#' Get names of the reads showing alternative allele of a mutation
#'
#' Extract the names of the reads in a bam file that support the variant allele of a single mutation
#' @param bam path to bam file
#' @param chr Chromosome name for the mutation
#' @param pos Chromosome position for the mutation
#' @param alt The alternative allele of the mutation
#' @param tag the RG tag if the bam has more than one samplee
#' @return A character vector having the read names
#' @export 

get_mutation_read_names <- function(bam, chr, pos, alt, tag = "") {

    assertthat::assert_that(!missing(bam), is.character(bam), length(bam) == 1, file.exists(bam))

    assertthat::assert_that(is.character(tag), length(tag) == 1)

    assertthat::assert_that(is.numeric(pos), pos > 0, pos %% 1 == 0, length(pos) == 1)

    assertthat::assert_that(is.character(chr), length(chr) == 1, chr %in% get_bam_chr(bam))

    assertthat::assert_that(is.character(alt), length(alt) == 1, alt %in% c("C","T","G","A"))

    gr <- GenomicRanges::GRanges(chr, IRanges::IRanges(pos, pos))
  
    if(tag == ""){
   
      stackedStrings <- GenomicAlignments::stackStringsFromBam(bam, use.names = T, param = gr)
  
    } else {
  
  	  assertthat::assert_that(verify_tag(bam = bam, tag = tag), 
          msg = "Specified tag not found")

      stackedStrings <- GenomicAlignments::stackStringsFromBam(bam, use.names = T, 
          param = Rsamtools::ScanBamParam(tagFilter = list("RG" = tag), which = gr)) 
  
    }

    sm <- get_bam_SM(bam = bam, tag = tag)
    if(length(stackedStrings) != 0){
      
      out <- data.frame(ID = paste(sm, names(stackedStrings), sep = "_"), 
  	    seq = as.data.frame(stackedStrings)[,1])
  
      return(unique(as.character(out[out$seq == alt, "ID"])))
   
   } else {

      return(character(0))
   
   }
}


#' Get names of the reads showing alternative allele of a mutation
#'
#' Extract the names of the reads in a bam file that support the variant allele of a single mutation
#' @param bam path to bam file
#' @param mutations A data frame containing the mutations. Must have the columns CHROM, POS, REF, ALT.
#' @param tag the RG tag if the bam has more than one samplee
#' @return A list with length equal to the number of mutations. Each element is a character vector with the read names.
#' @export 

get_mutations_read_names <- function(bam, mutations, tag = "") {
  
    assertthat::assert_that(!missing(bam), is.character(bam), 
      length(bam) == 1, file.exists(bam))

    assertthat::assert_that(!missing(mutations), is.data.frame(mutations), assertthat::not_empty(mutations), 
        assertthat::has_name(mutations, c("CHROM", "POS", "REF", "ALT")))

    assertthat::assert_that(all(nchar(mutations$REF) == 1),
        all(nchar(mutations$ALT) == 1),
        msg = "Only SNVs are supported") 
    
    assertthat::assert_that(is.character(mutations$REF), is.character(mutations$ALT),
        all(mutations$REF %in% c("A", "C", "T", "G")), all(mutations$ALT %in% c("A", "C", "T", "G")),
        msg = "REF and ALT in mutations should be characters having basepairs")

    assertthat::assert_that(is.numeric(mutations$POS), all(mutations$POS > 0))

    assertthat::assert_that(all(mutations$CHROM %in% get_bam_chr(bam)))

    read_names <- purrr::pmap(list(mutations$CHROM, mutations$POS, mutations$ALT),
        function(chr, pos, alt) {
            counts <- get_mutation_read_names(bam = bam, tag = tag,
              chr = chr, pos = pos, alt = alt)
        })
    
    names(read_names) <- paste0(mutations$CHROM, ":", mutations$POS, "_", 
        mutations$REF, "_", mutations$ALT)
    
    return(read_names)

}