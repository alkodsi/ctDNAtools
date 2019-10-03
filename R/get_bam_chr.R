#' Helper function to extract chromosome names from a bam file

get_bam_chr <- function(bam) {
   
  header <- Rsamtools::scanBamHeader(bam)[[1]]$text
  chr <- gsub("SN:", "", purrr::map_chr(header[ names(header) == "@SQ" ], 1))
  return(chr)
   
}