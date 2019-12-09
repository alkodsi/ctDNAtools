#' The function summarizes fragment size in defined genomic regions
#'
#' @param bam the input bam file
#' @param regions data frame contatining the genomic regions. Must have the columns chr, start and end.
#' @param tag  the RG tag if the bam has more than one samplee.
#' @param isProperPair a logical wheter to return only proper pairs (T), only improper pairs (F), or it does not matter (NA).
#' @param mapqFilter mapping quality threshold for considering reads.
#' @param min_size Integer with the lowest fragment length.
#' @param max_size Integer with the highest fragment length.
#' @param ignore_trimmed logical, whether to remove reads that have been hard trimmed.
#' @param different_strands logical, whether to keep only reads whose mates map to different strand.
#' @param simple_cigar logical, whether to include only reads with simple cigar.
#' @param summary_functions a named list containing the R functions used for summarization, e.g. mean, sd.
#' @return a data frame with the first column having the regions in the format of chr:start-end, and other columns correspond to summary_functions.
#' @export
#' @importFrom rlang .data
#' @importFrom magrittr %>%

summarize_fragment_size <- function(bam, regions, tag = "", isProperPair = NA, mapqFilter = 30, 
    different_strands = T, min_size = 1, max_size = 400, simple_cigar = F,  ignore_trimmed = T,
    summary_functions = list(Mean = mean, Median = median)) {
   
    assertthat::assert_that(assertthat::has_name(regions, c("chr", "start", "end")),
        is.data.frame(regions))

    assertthat::assert_that(!missing(bam), is.character(bam), length(bam) == 1, file.exists(bam))

    assertthat::assert_that(is.character(tag), length(tag) == 1)

    assertthat::assert_that(is.logical(isProperPair), is.logical(ignore_trimmed), 
        is.logical(different_strands), is.logical(simple_cigar),
        is.numeric(min_size), is.numeric(max_size), max_size > min_size, 
        length(isProperPair) == 1, length(ignore_trimmed) == 1, length(different_strands) == 1, 
        length(simple_cigar) == 1, length(min_size) == 1, length(max_size) == 1)

    assertthat::assert_that(is.list(summary_functions), !is.null(names(summary_functions)),
        msg = "summary_functions must be a named list")

    assertthat::assert_that(all(map_lgl(list(mean = mean, median = median), is.function)),
        msg = "all elements of summary_functions must be functions")

    assertthat::assert_that(all(purrr::map_dbl(purrr::invoke_map(summary_functions, x = c(1:5)),length)==1),
        msg = "Functions in summary_functions must produce output of length 1")

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

    regions$ID <- paste0(regions$chr, ":", regions$start, "-", regions$end)

    regions_gr <- GenomicRanges::GRanges(regions$chr, IRanges::IRanges(regions$start, regions$end))

    overlaps <- as.data.frame(GenomicRanges::findOverlaps(reads, regions_gr))

    fragments <- as.data.frame(reads) %>%
        dplyr::mutate(Region = NA, size = abs(.data$isize.first)) %>%
        dplyr::select(Fragment = qname.first, size, Region)

    fragments[overlaps[,1], "Region"] <- regions[overlaps[,2], "ID"]
   
    summary <- fragments %>% 
        tidyr::replace_na(list(Region = "offTargets")) %>%
        dplyr::group_by(.data$Region) %>%
        dplyr::summarize_if(is.numeric, summary_functions) %>%
        as.data.frame()

    colnames(summary)[-1] <- paste(sm, names(summary_functions), sep = "_")

    return(summary) 
}