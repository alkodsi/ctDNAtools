#' Get names of the reads showing alternative allele of a mutation
#'
#' Extract the names of the reads in a bam file that support the variant allele of a single mutation
#' @param bam path to bam file
#' @param mutations A data frame containing the mutations. Must have the columns CHROM, POS, REF, ALT.
#' @param what Either Alt: names of the reads that exhibit variant alleles in mutations, 
#' or Ref: names of the reads the exhibit reference alleles in mutations
#' @param tag the RG tag if the bam has more than one samplee
#' @return A list with length equal to the number of mutations. Each element is a character vector with the read names.
#' @export 

get_mutations_read_names <- function(bam, mutations, what = c("ALT","REF"), min_base_quality = 20, tag = "") {
    
    assertthat::assert_that(!missing(bam), is.character(bam),
     length(bam) == 1, file.exists(bam))
    
    assertthat::assert_that(!missing(mutations), 
        is.data.frame(mutations), assertthat::not_empty(mutations), 
        assertthat::has_name(mutations, c("CHROM", "POS", "REF", "ALT")))
    
    assertthat::assert_that(all(nchar(mutations$REF) == 1), 
        all(nchar(mutations$ALT) == 1), msg = "Only SNVs are supported")
    
    assertthat::assert_that(is.character(mutations$REF), is.character(mutations$ALT), 
        all(mutations$REF %in% c("A", "C", "T", "G")), 
        all(mutations$ALT %in% c("A", "C", "T", "G")),
        msg = "REF and ALT in mutations should be characters having basepairs")
    
    assertthat::assert_that(is.numeric(mutations$POS), all(mutations$POS > 0))
    
    assertthat::assert_that(all(mutations$CHROM %in% get_bam_chr(bam)))
    
    what <- match.arg(what)

    read_names <- purrr::pmap(list(mutations$CHROM, mutations$POS, mutations[,what]), 
        function(chr, pos, alt) {
            counts <- get_mutation_read_names(bam = bam, tag = tag, chr = chr, pos = pos, 
                alt = alt, min_base_quality = min_base_quality)
        })
    
    names(read_names) <- paste0(mutations$CHROM, ":", mutations$POS, "_", mutations$REF, "_", mutations$ALT)
    
    return(read_names)
    
}



