
#' @export
create_black_list <- function(bam_list, targets, reference, vaf_threshold = 0.1, bam_list_tags = rep("",length(bam_list)), 
    min_base_quality = 20, max_depth = 1e+05, min_mapq = 30, mean_vaf_quantile = 0.95, 
    min_samples_one_read = max(2, ceiling(length(bam_list) * 0.3)), 
    min_samples_two_reads = max(2, ceiling(length(bam_list) * 0.05))) {
    
    assertthat::assert_that(class(reference) == "BSgenome")
    assertthat::assert_that(is.data.frame(targets), assertthat::not_empty(targets), 
        assertthat::has_name(targets, c("chr", "start", "end")))
        
    assertthat::assert_that(is.character(bam_list), all(file.exists(bam_list)))
        
    bam_list_chr <- purrr::map(bam_list, get_bam_chr)
        
    assertthat::assert_that(all(purrr::map_lgl(bam_list_chr, ~ all(.x %in% GenomeInfoDb::seqnames(reference)))), 
        msg = "The chromosomes in at least one of the specified bams in bam_list don't match the reference")
        
    assertthat::assert_that(is.character(bam_list_tags), 
        length(bam_list_tags) == length(bam_list))
        
    if (any(bam_list_tags != "")) {
            
        tag_verification <- purrr::map2_lgl(bam_list[bam_list_tags != ""],
            bam_list_tags[bam_list_tags != ""], verify_tag)
            
        assertthat::assert_that(all(tag_verification), 
            msg = paste("specified tag for", paste(bam_list[bam_list_tags != ""][!tag_verification], collapse = " , "), 
           "is not correct"))    
    }

    background_panel <- purrr::map2(bam_list, bam_list_tags,
   	  ~ create_black_list_instance(bam = .x, tag = .y, targets = targets, reference = reference,
   	  	   vaf_threshold = vaf_threshold, min_base_quality = min_base_quality, min_mapq = min_mapq,
   	  	   max_depth = max_depth))
    
    sm <- make.unique(purrr::map_chr(bam_list, get_bam_SM))

    depth <- purrr::reduce(purrr::map(background_panel, "depth"), dplyr::full_join, by = "Locus") %>%
        rlang::set_names(c("Locus", sm))

    alt <- purrr::reduce(purrr::map(background_panel, "alt"), dplyr::full_join, by = "Locus") %>%
        rlang::set_names(c("Locus", sm))
    
    vaf <- purrr::reduce(purrr::map(background_panel, "vaf"), dplyr::full_join, by = "Locus") %>%
        rlang::set_names(c("Locus", sm))
    
    mean_vaf <- rowMeans(vaf[, -1, drop = F], na.rm = T)
    mean_vaf_idx <- which(mean_vaf >= quantile(mean_vaf, mean_vaf_quantile, na.rm = T))

    samples_one_read <- rowSums(alt[, -1, drop = F] >= 1, na.rm =T)
    samples_one_read_idx <- which(samples_one_read >= min_samples_one_read)

    samples_two_reads <- rowSums(alt[, -1, drop = F] >=2, na.rm = T)
    samples_two_reads_idx <- which(samples_two_reads >= min_samples_two_reads)
    
    idx <- unique(c(mean_vaf_idx, samples_one_read_idx, samples_two_reads_idx))

    black_list <- depth$Locus[idx]

    return(list(black_list = black_list, depth = depth, alt = alt, vaf = vaf))
}


create_black_list_instance <- function(bam, targets, reference, vaf_threshold = 0.1, tag = "", 
    min_base_quality = 20, max_depth = 1e+05, min_mapq = 30) {

    gr <- GenomicRanges::reduce(GenomicRanges::GRanges(targets$chr, IRanges::IRanges(targets$start, targets$end)))

    if (tag == "") {
    
        sbp <- Rsamtools::ScanBamParam(which = gr)
    
    } else {
    
        sbp <- Rsamtools::ScanBamParam(which = gr, tagFilter = list(RG = tag))
    
    }

    pileupParam <- Rsamtools::PileupParam(max_depth = max_depth, min_base_quality = min_base_quality, 
        min_mapq = min_mapq, distinguish_strands = F, include_deletions = F, include_insertions = F)

    p <- Rsamtools::pileup(bam, scanBamParam = sbp, pileupParam = pileupParam) %>% 
        tidyr::pivot_wider(names_from = .data$nucleotide, values_from = .data$count, 
          values_fill = list(count = 0)) %>% 
        as.data.frame() %>% 
        dplyr::mutate(ref = as.character(BSgenome::getSeq(reference, 
            GenomicRanges::GRanges(.data$seqnames, IRanges::IRanges(.data$pos, .data$pos))))) %>% 
        dplyr::mutate(depth = .data$A + .data$C + .data$G + .data$T)
    
    pAnn <- dplyr::mutate(p, refCount = purrr::map2_dbl(c(1:nrow(p)), p$ref, ~p[.x, .y]),
        nonRefCount = .data$depth - .data$refCount, vaf = .data$nonRefCount/.data$depth) %>% 
        dplyr::mutate(nonRefCount = ifelse(vaf >= vaf_threshold, NA, nonRefCount),
        	depth = ifelse(vaf >= vaf_threshold, NA, depth),
            vaf = ifelse(vaf >= vaf_threshold, NA, vaf), 
            Locus = paste(.data$seqnames, .data$pos, sep = "_")) %>%
        dplyr::select(.data$Locus, .data$depth, .data$nonRefCount, .data$vaf)

    out <- list(depth = data.frame(Locus = pAnn$Locus, depth = pAnn$depth, stringsAsFactors = F),
                alt = data.frame(Locus = pAnn$Locus, alt = pAnn$nonRefCount, stringsAsFactors = F),
                vaf = data.frame(Locus = pAnn$Locus, vaf = pAnn$vaf, stringsAsFactors = F))
    return(out)
}