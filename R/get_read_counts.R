#' Get read counts for a specific base in the genome
#'
#' Uses samtools pileup to get the read counts for each base in the genomic position specified
#' @param chr chromosome name
#' @param pos genomic coordinate
#' @param bam path to bam file
#' @param tag the RG tag if the bam has more than one sample
#' @param min_base_quality minimum base quality for a read to be counted
#' @param max_depth maximum depth above which sampling will happen
#' @param min_mapq the minimum mapping quality for a read to be counted
#' @return a list, number of reads for each of the four basepairs

get_read_counts <- function(chr, pos, bam, tag = "", min_base_quality = 20,
     max_depth = 1e+05,  min_mapq = 30) {

    gr <- GenomicRanges::GRanges(chr, IRanges::IRanges(pos, pos))
    
    if (tag == "") {
        sbp <- Rsamtools::ScanBamParam(which = gr)
    } else {
        sbp <- Rsamtools::ScanBamParam(which = gr, tagFilter = list(RG = tag))
    }

    pileupParam <- Rsamtools::PileupParam(max_depth = max_depth, 
        min_base_quality = min_base_quality, min_mapq = min_mapq, 
        distinguish_strands = F, include_deletions = F, include_insertions = F)
    
    p <- Rsamtools::pileup(bam, scanBamParam = sbp, pileupParam = pileupParam)

    cbase <- ifelse("C" %in% p$nucleotide, p[p$nucleotide == "C", "count"], 0)
    gbase <- ifelse("G" %in% p$nucleotide, p[p$nucleotide == "G", "count"], 0)
    abase <- ifelse("A" %in% p$nucleotide, p[p$nucleotide == "A", "count"], 0)
    tbase <- ifelse("T" %in% p$nucleotide, p[p$nucleotide == "T", "count"], 0)
   
    return(list(A = abase, C = cbase, G = gbase, T = tbase))
}

#' Helper function to verify a tag exists
verify_tag <- function(bam, tag) {

   header <- Rsamtools::scanBamHeader(bam) 
   header_text <- header[[1]]$text
   rg <- header_text[ names(header_text) == "@RG" ]
   tags <- gsub("ID:", "", purrr::map_chr(rg, 1))
   return(tag %in% tags)

}

#' Helper function to extract SM from a bam file
get_bam_SM <- function(bam, tag = "") {

    header <- Rsamtools::scanBamHeader(bam)
    if(tag == ""){
        rg <- header[[1]]$text$`@RG`
        sm <- gsub("SM:", "", rg[ grepl("SM", rg) ])
    } else {
        header_text <- header[[1]]$text
        rg <- header_text[ names(header_text) == "@RG" ]
        rg <- rg[[ which(grepl(tag, purrr::map_chr(rg, 1))) ]]
        sm <- gsub("SM:", "", rg[ grepl("SM",rg) ])
    }
   return(sm) 

}

#' Helper function to extract chromosome names from a bam file
get_bam_chr <- function(bam) {
   
  header <- Rsamtools::scanBamHeader(bam)[[1]]$text
  chr <- gsub("SN:", "", purrr::map_chr(header[ names(header) == "@SQ" ], 1))
  return(chr)
   
}