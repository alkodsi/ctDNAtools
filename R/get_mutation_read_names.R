
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
    
    assertthat::assert_that(is.numeric(pos), pos > 0, pos%%1 == 0, length(pos) == 1)
    
    assertthat::assert_that(is.character(chr), length(chr) == 1, chr %in% get_bam_chr(bam))
    
    assertthat::assert_that(is.character(alt), length(alt) == 1, alt %in% c("C", "T", "G", "A"))
    
    gr <- GenomicRanges::GRanges(chr, IRanges::IRanges(pos, pos))
    
    if (tag == "") {
        
        stackedStrings <- GenomicAlignments::stackStringsFromBam(bam, use.names = T, 
            param = gr)
        
    } else {
        
        assertthat::assert_that(verify_tag(bam = bam, tag = tag), msg = "Specified tag not found")
        
        stackedStrings <- GenomicAlignments::stackStringsFromBam(bam, use.names = T, 
            param = Rsamtools::ScanBamParam(tagFilter = list(RG = tag), which = gr))
        
    }
    
    sm <- get_bam_SM(bam = bam, tag = tag)
    if (length(stackedStrings) != 0) {
        
        out <- data.frame(ID = paste(sm, names(stackedStrings), sep = "_"), 
            seq = as.data.frame(stackedStrings)[, 1],
            stringsAsFactors = F)
        
        return(unique(out[out$seq == alt, "ID"]))
        
    } else {
        
        return(character(0))
        
    }
}
