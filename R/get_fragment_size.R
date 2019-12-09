#' Get fragment lengths from bam
#'
#' A function to extract fragment lengths from a bam file. Optionally, given a mutation data frame, it can categorize read lengths 
#' in mutated vs non-mutated reads.
#' @param bam path to bam file.
#' @param mutations An optional data frame with mutations. Must have the columns CHROM, POS, REF, ALT.
#' @param tag the RG tag if the bam has more than one samplee.
#' @param isProperPair a logical wheter to return only proper pairs (T), only improper pairs (F), or it does not matter (NA).
#' @param mapqFilter mapping quality threshold for considering reads.
#' @param min_size Integer with the lowest fragment length.
#' @param max_size Integer with the highest fragment length.
#' @param ignore_trimmed logical, whether to remove reads that have been hard trimmed.
#' @param different_strands logical, whether to keep only reads whose mates map to different strand.
#' @param simple_cigar logical, whether to include only reads with simple cigar.
#' @return A data frame with the columns Sample (SM tag in bam), ID (read ID), size (fragment size), and category (only if mutations is provided).
#' @export 
#' @importFrom rlang .data
#' @importFrom magrittr %>%

get_fragment_size <- function(bam, mutations = NULL, tag = "", isProperPair = NA, mapqFilter = 30,
	min_size = 1, max_size = 400, ignore_trimmed = T, different_strands = T, simple_cigar = F) {

    assertthat::assert_that(!missing(bam), is.character(bam), length(bam) == 1, file.exists(bam))

    assertthat::assert_that(is.character(tag), length(tag) == 1)
    
    assertthat::assert_that(is.logical(isProperPair), is.logical(ignore_trimmed), 
    	is.logical(different_strands), is.logical(simple_cigar),
    	is.numeric(min_size), is.numeric(max_size), max_size > min_size, 
    	length(isProperPair) == 1, length(ignore_trimmed) == 1, length(different_strands) == 1, 
    	length(simple_cigar) == 1, length(min_size) == 1, length(max_size) == 1)

    if(!is.null(mutations)) {
        
        assertthat::assert_that(is.data.frame(mutations), assertthat::not_empty(mutations), 
        
        assertthat::has_name(mutations, c("CHROM", "POS", "REF", "ALT")))

        assertthat::assert_that(all(nchar(mutations$REF) == 1),
            all(nchar(mutations$ALT) == 1),
            msg = "Only SNVs are supported") 
    
        assertthat::assert_that(is.character(mutations$REF), is.character(mutations$ALT),
            all(mutations$REF %in% c("A", "C", "T", "G")), all(mutations$ALT %in% c("A", "C", "T", "G")),
            msg = "REF and ALT in mutations should be characters having basepairs")

        assertthat::assert_that(is.numeric(mutations$POS), all(mutations$POS > 0))

        assertthat::assert_that(all(mutations$CHROM %in% get_bam_chr(bam)))
    }

    flag <- Rsamtools::scanBamFlag(isPaired = T, isProperPair = isProperPair, isUnmappedQuery = F,
              hasUnmappedMate = F, isSecondaryAlignment = F,
              isNotPassingQualityControls = F, isDuplicate = F, isSupplementaryAlignment = F)
    
    if(tag != ""){
	
	    assertthat::assert_that(verify_tag(bam = bam, tag = tag), 
           msg = "Specified tag not found")
	    
	    sbp <- Rsamtools::ScanBamParam(flag = flag, mapqFilter = mapqFilter, simpleCigar = simple_cigar,
	       tagFilter = list(RG = tag), what = c("qname","flag","qwidth","isize"))
    
    } else {

        sbp <- Rsamtools::ScanBamParam(flag = flag, mapqFilter = mapqFilter, simpleCigar = simple_cigar,
           what = c("qname","flag","qwidth","isize"))
    
    }
    
    sm <- get_bam_SM(bam = bam, tag = tag)
    
    reads <- GenomicAlignments::readGAlignmentPairs(file = bam, param = sbp)
    
    reads <- reads[abs(reads@first@elementMetadata$isize) > min_size & 
        abs(reads@first@elementMetadata$isize) < max_size]

    if(different_strands) {

        reads <- reads[xor(BiocGenerics::strand(reads@first) == "+",
            BiocGenerics::strand(reads@last) == "+")]
    
    }

    if(ignore_trimmed){

        max_qwidth <- max(reads@first@elementMetadata$qwidth)

        reads <- reads[reads@first@elementMetadata$qwidth == max_qwidth & 
            reads@last@elementMetadata$qwidth == max_qwidth]
   
    }

    fragment_lengths <- as.data.frame(reads) %>%
        dplyr::mutate(Sample = sm,
            ID = paste(sm, .data$qname.first, sep = "_"),
            chr = .data$seqnames.first,
            start = pmin(.data$start.first, .data$end.first, .data$start.last, .data$end.last),
            end = pmax(.data$start.first, .data$end.first, .data$start.last, .data$end.last),
            size = abs(.data$isize.first),
            stringsAsFactors = F) %>%
        dplyr::select(.data$Sample, .data$ID, .data$chr,
            .data$start, .data$end, .data$size)

    if(!is.null(mutations)) {
     
        read_names <- unique(unlist(purrr::map(get_mutations_read_names(bam = bam, tag = tag, mutations = mutations), "alt")))
        
        fragment_lengths <- dplyr::mutate(fragment_lengths, 
        	category = ifelse(.data$ID %in% read_names, "mutated", "other"))
    }

    return(fragment_lengths)
}

#' @export
plot_fragment_size_region <- function(fragments1, fragments2 = NULL, chr, start, end, 
    smoothing_window = 100, genome = c("hg19","hg38")) {
    
    fragments1_bychr <- split(fragments1, fragments1$chr)
    
    fragments1 <- purrr::map_dfr(fragments1_bychr, 
        ~ dplyr::mutate(.x, smoothed_size = RcppRoll::roll_mean(.data$size, n = smoothing_window, fill = 0))) %>%
       dplyr::filter(!is.na(.data$smoothed_size))
    
    sample_levels <- unique(fragments1$Sample)

    if(!is.null(fragments2)){
        fragments2_bychr <- split(fragments2, fragments2$chr)
         
        fragments2<- purrr::map_dfr(fragments2_bychr, 
             ~ dplyr::mutate(.x, smoothed_size = RcppRoll::roll_mean(.data$size, n = smoothing_window, fill = 0))) %>%
               dplyr::filter(!is.na(.data$smoothed_size))
        
        sample_levels <- c(sample_levels, unique(fragments2$Sample))
    
    }
    
    genome <- match.arg(genome)

    gtrack <- Gviz::GenomeAxisTrack()

    itrack <- Gviz::IdeogramTrack( chromosome=chr, genome = genome)
    

    data1 <- GenomicRanges::GRanges(fragments1$chr, IRanges(fragments1$pos, fragments1$pos), size = fragments1$smoothed_size)

    covTrack1 <- Gviz::DataTrack(as(coverage(data1), "GRanges"), genome = genome,
         type = "l", name = "Coverage", groups=factor(unique(fragments1$Sample), levels=sample_levels))

    dtrack1 <- Gviz::DataTrack(data1,  genome = genome, name = "Fragment Length",
         groups=factor(unique(fragments1$Sample), levels=sample_levels))
    
    if(!is.null(fragments2)){
   
        data2 <- GenomicRanges::GRanges(fragments2$chr, IRanges::IRanges(fragments2$pos, fragments2$pos), size2 = fragments2$smoothed_size)
        
        covTrack2 <- DataTrack(as(coverage(data2), "GRanges"), genome = genome,
           type = "l", name = "Coverage", groups=factor(unique(fragments2$Sample), levels=sample_levels))
        
        dtrack2 <- Gviz::DataTrack(data2, groups=factor(unique(fragments2$Sample), levels=sample_levels),
             genome = genome, name = "Fragment Length")

        overlay <- Gviz::OverlayTrack(trackList = list(dtrack1, dtrack2))
   
        overlay2 <- Gviz::OverlayTrack(trackList = list(covTrack1, covTrack2))
        #ylims <- extendrange(range(c(values(dtrack1), values(dtrack2))))
   
        Gviz::plotTracks(list(itrack, gtrack, overlay2, overlay), from = start, to = end, chromosome = chr,  type="a",
            background.panel="#FFFEDB", background.title="darkblue")
   
    } else {

        Gviz::plotTracks(list(itrack, gtrack, covTrack1, dtrack1), from = start, to = end, chromosome = chr, type="a",
            background.panel="#FFFEDB", background.title="darkblue")

    }

}