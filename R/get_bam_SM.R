#' Helper function to extract SM from a bam file
#' @param bam Path to the bam file
#' @param tag The RG tag for the bam file
#' @keywords internal

get_bam_SM <- function(bam, tag = "") {
  header <- Rsamtools::scanBamHeader(bam)
  if (tag == "") {
    rg <- header[[1]]$text$`@RG`
    sm <- gsub("SM:", "", rg[grepl("SM", rg)])
  } else {
    header_text <- header[[1]]$text
    rg <- header_text[names(header_text) == "@RG"]
    rg <- rg[[which(grepl(tag, purrr::map_chr(rg, 1)))]]
    sm <- gsub("SM:", "", rg[grepl("SM", rg)])
  }
  if (is.null(rg)) {
    sm <- basename(bam)
  }

  return(sm)
}
