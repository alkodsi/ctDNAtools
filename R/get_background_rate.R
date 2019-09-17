#' Compute the background mismatch rate from a bam file
#'
#' Runs through the target regions base by base counting the mismatches. Then it divides sum(msimatches)/sum(depths) for all bases in the targets
#' @param bam path to bam file
#' @param targets a data frame with the target regions. Must have three columns: chr, start and end
#' @param reference the reference genome in BSgenome format
#' @param vaf_threshold the bases with higher than this VAF threshold will be ignored in the calculation (real mutations)
#' @param tag the RG tag if the bam has more than one sample
#' @param min_base_quality minimum base quality for a read to be counted
#' @param max_depth maximum depth above which sampling will happen
#' @param min_mapq the minimum mapping quality for a read to be counted
#' @return a list containing the general mismatch rate and substitution-specific rates
#' @export

get_background_rate <- function(bam, targets, reference, vaf_threshold = 0.1, tag = "", 
    min_base_quality = 20, max_depth = 1e+05, min_mapq = 30) {
    

    gr <- GRanges(targets$chr, IRanges(targets$start, targets$end))
    
    assertthat::assert_that(class(reference) == "BSgenome")
    assertthat::assert_that(is.data.frame(targets), 
          assertthat::not_empty(targets), 
          assertthat::has_name(targets, c("chr", "start", "end")))
    
    if (tag == "") {
        sbp <- Rsamtools::ScanBamParam(which = gr)
    } else {
        sbp <- Rsamtools::ScanBamParam(which = gr, tagFilter = list(RG = tag))
    }

    pileupParam <- Rsamtools::PileupParam(max_depth = max_depth, min_base_quality = min_base_quality, 
        min_mapq = min_mapq, distinguish_strands = F, include_deletions = F, 
        include_insertions = F)
    
    p <- Rsamtools::pileup(bam, scanBamParam = sbp, pileupParam = pileupParam) %>% 
        tidyr::spread(-nucleotide, count, fill = 0) %>% 
        dplyr::mutate(ref = as.character(getSeq(reference, GRanges(seqnames, IRanges(pos, pos))))) %>% 
        dplyr::mutate(depth = A + C + G + T)
    
    pAnn <- dplyr::mutate(p, refCount = purrr::map2_dbl(c(1:nrow(p)), p$ref, ~p[.x, .y]),
                      nonRefCount = depth - refCount) %>% 
            dplyr::filter(nonRefCount/depth < vaf_threshold)
    
    rate_by_sub <- pAnn %>% 
          dplyr::group_by(ref) %>% 
          dplyr::summarize(
              depth = sum(as.numeric(depth)), 
              A = sum(as.numeric(A)), 
              C = sum(as.numeric(C)), 
              G = sum(as.numeric(G)), 
              T = sum(as.numeric(T))) %>%
           as.data.frame()
    
    CA <- (rate_by_sub[ rate_by_sub$ref == "C", "A"] + rate_by_sub[ rate_by_sub$ref == "G", "T"]) /
          (rate_by_sub[ rate_by_sub$ref == "C", "depth"] + rate_by_sub[ rate_by_sub$ref == "G", "depth"])
    
    CG <- (rate_by_sub[ rate_by_sub$ref == "C", "G"] + rate_by_sub[ rate_by_sub$ref == "G", "C"]) /
          (rate_by_sub[ rate_by_sub$ref == "C", "depth"] + rate_by_sub[ rate_by_sub$ref == "G", "depth"])              
    
    CT <- (rate_by_sub[ rate_by_sub$ref == "C", "T"] + rate_by_sub[ rate_by_sub$ref == "G", "A"]) /
          (rate_by_sub[ rate_by_sub$ref == "C", "depth"] + rate_by_sub[ rate_by_sub$ref == "G", "depth"])

    TA <- (rate_by_sub[ rate_by_sub$ref == "T", "A"] + rate_by_sub[ rate_by_sub$ref == "A", "T"]) /
          (rate_by_sub[ rate_by_sub$ref == "T", "depth"] + rate_by_sub[ rate_by_sub$ref == "A", "depth"])
    
    TC <- (rate_by_sub[ rate_by_sub$ref == "T", "C"] + rate_by_sub[ rate_by_sub$ref == "A", "G"]) /
          (rate_by_sub[ rate_by_sub$ref == "T", "depth"] + rate_by_sub[ rate_by_sub$ref == "A", "depth"])              
    
    TG <- (rate_by_sub[ rate_by_sub$ref == "T", "G"] + rate_by_sub[ rate_by_sub$ref == "A", "C"]) /
          (rate_by_sub[ rate_by_sub$ref == "T", "depth"] + rate_by_sub[ rate_by_sub$ref == "A", "depth"])
          
    rate <- sum(as.numeric(pAnn$nonRefCount))/sum(as.numeric(pAnn$depth))
    
    return(list(rate = rate, CA = CA, CG = CG, CT = CT,
                TA = TA, TC = TC, TG = TG))
}

