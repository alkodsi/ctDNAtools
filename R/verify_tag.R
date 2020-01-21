#' Helper function to verify a tag exists

#' @param bam path to bam file
#' @param tag the RG tag if the bam has more than one sample
#' @return logical, whether tag in bam
#' @keywords internal

verify_tag <- function(bam, tag) {
  header <- Rsamtools::scanBamHeader(bam)
  header_text <- header[[1]]$text
  rg <- header_text[names(header_text) == "@RG"]
  tags <- gsub("ID:", "", purrr::map_chr(rg, 1))
  return(tag %in% tags)
}
