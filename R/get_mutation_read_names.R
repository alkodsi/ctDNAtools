
#' Get names of the reads showing alternative allele of a mutation
#'
#' Extract the names of the reads in a bam file that support the variant allele of a single mutation
#' @param bam path to bam file
#' @param chr Chromosome name for the mutation
#' @param pos Chromosome position for the mutation
#' @param alt The alternative allele of the mutation
#' @param tag the RG tag if the bam has more than one samplee
#' @return A character vector having the read names
#' @importFrom rlang .data
#' @export 

get_mutation_read_names <- function(bam, chr, pos, alt, tag = "", min_base_quality = 20) {
    
    assertthat::assert_that(!missing(bam), is.character(bam), length(bam) == 1, file.exists(bam))
    
    assertthat::assert_that(is.character(tag), length(tag) == 1)
    
    assertthat::assert_that(is.numeric(pos), pos > 0, pos%%1 == 0, length(pos) == 1)
    
    assertthat::assert_that(is.character(chr), length(chr) == 1, chr %in% get_bam_chr(bam))
    
    assertthat::assert_that(is.character(alt), length(alt) == 1, alt %in% c("C", "T", "G", "A"))
    
    gr <- GenomicRanges::GRanges(chr, IRanges::IRanges(pos, pos))
    
    if (tag == "") {
        
        sbp <- Rsamtools::ScanBamParam(flag = Rsamtools::scanBamFlag(isUnmappedQuery = F, 
            isSecondaryAlignment = F, isNotPassingQualityControls = F, isDuplicate = F),
            which = gr)

        stackedStrings <- GenomicAlignments::stackStringsFromBam(bam, use.names = T, 
            param = sbp)

         stackedQuals <- GenomicAlignments::stackStringsFromBam(bam, use.names = T,
             what = "qual", param = sbp)
        
    } else {
        
        assertthat::assert_that(verify_tag(bam = bam, tag = tag), msg = "Specified tag not found")
        
        sbp <- Rsamtools::ScanBamParam(flag = Rsamtools::scanBamFlag(isUnmappedQuery = F, 
            isSecondaryAlignment = F, isNotPassingQualityControls = F, isDuplicate = F),
            which = gr, tagFilter = list(RG = tag))

        stackedStrings <- GenomicAlignments::stackStringsFromBam(bam, use.names = T, 
            param = sbp)
        
        stackedQuals <- GenomicAlignments::stackStringsFromBam(bam, use.names = T, 
            what = "qual", param = sbp)
    }
    
    sm <- get_bam_SM(bam = bam, tag = tag)

    if (length(stackedStrings) != 0) {
        
        out <- data.frame(ID = paste(sm, names(stackedStrings), sep = "_"), 
            seq = as.data.frame(stackedStrings)[, 1],
            qual = as.data.frame(methods::as(Biostrings::PhredQuality(stackedQuals), "IntegerList"))$value,
            stringsAsFactors = F) %>%
            dplyr::filter(.data$qual >= min_base_quality)
        
        return(unique(out[out$seq == alt, "ID"]))
        
    } else {
        
        return(character(0))
        
    }
}
