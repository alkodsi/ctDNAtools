#' Helper function to extract chromosome names from a bam file
#' @param bam the path to bam file
#' @keywords internal

get_bam_chr <- function(bam) {
  chr <- GenomeInfoDb::seqnames(GenomeInfoDb::seqinfo(Rsamtools::BamFile(bam)))

  return(chr)
}
