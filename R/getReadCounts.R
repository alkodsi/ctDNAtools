#' Get read counts for a specific base in the genome
#'
#' Uses samtools pileup to get the read counts of the genomic position and the base specified
#' @param chr chromosome name
#' @param pos genomic coordinate
#' @param base the basepair for which count the reads
#' @param bam path to bam file
#' @param tag the RG tag if the bam has more than one sample
#' @param min_base_quality minimum base quality for a read to be counted
#' @param max_depth maximum depth above which sampling will happen
#' @param include_indels whether to include indels in the pileup
#' @param min_mapq the minimum mapping quality for a read to be counted
#' @return a scalar, number of reads matching the base in the specified genomic position

getReadCounts <- function(chr, pos, base, bam, tag = "", min_base_quality = 20, max_depth = 1e+05, 
    include_indels = F, min_mapq = 30) {

    gr <- GRanges(chr, IRanges(pos, pos))
    
    if (tag == "") {
        sbp <- Rsamtools::ScanBamParam(which = gr)
    } else {
        sbp <- Rsamtools::ScanBamParam(which = gr, tagFilter = list(RG = tag))
    }

    pileupParam <- Rsamtools::PileupParam(max_depth = max_depth, min_base_quality = min_base_quality, 
        min_mapq = min_mapq, distinguish_strands = F, include_deletions = include_indels, 
        include_insertions = include_indels)
    
    p <- Rsamtools::pileup(bam, scanBamParam = sbp, pileupParam = pileupParam)

    depth <- ifelse(!base %in% p$nucleotide, 0, p[p$nucleotide == base, "count"])
    
    return(depth)
}
