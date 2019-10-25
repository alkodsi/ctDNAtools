
#' @export
#' @importFrom mclust mclustBIC

estimate_ctDNA_level <- function(mutations, cfDNA_conc = NA, bam, reference, targets, tag = "", vaf_column = "vaf", 
	ref_reads_column = "RefReads", alt_reads_column = "AltReads", background_rate = NA, remove_traces = T, 
	fdr_thr = 0.05, only_SNVs = F, use_clustering = T, override_selected_cluster = NA, ...) {
       
    assertthat::assert_that(is.character(ref_reads_column), is.character(alt_reads_column), 
    	is.character(vaf_column), length(ref_reads_column) == 1, length(alt_reads_column) == 1,
    	length(vaf_column) == 1)

    assertthat::assert_that(is.data.frame(mutations), assertthat::not_empty(mutations), 
        assertthat::has_name(mutations, c("CHROM", "POS", "REF", "ALT", vaf_column, ref_reads_column, alt_reads_column)))
    
    assertthat::assert_that(is.character(mutations$REF), is.character(mutations$ALT))
    
    assertthat::assert_that(is.numeric(mutations$POS), all(mutations$POS > 0))
    
    assertthat::assert_that(is.numeric(targets$start), all(mutations$start > 0), 
        is.numeric(targets$end), all(mutations$end > 0))
    
    assertthat::assert_that(is.logical(use_clustering), is.logical(only_SNVs), 
    	is.logical(remove_traces), length(use_clustering) == 1,  length(remove_traces) == 1,
    	length(only_SNVs) == 1, length(override_selected_cluster) == 1)
    
    mutations <- mutations[!is.na(mutations[,vaf_column]),]

    assertthat::assert_that(is.numeric(mutations[,vaf_column]), 
    	all(mutations[,vaf_column] <= 1), all(mutations[,vaf_column] >= 0))
    
    assertthat::assert_that(is.numeric(mutations[,ref_reads_column]), 
    	all(mutations[,ref_reads_column] %% 1 == 0), all(mutations[,ref_reads_column] %% 1 == 0),
    	msg = "specified ref_reads_column should have all integer numbers")

    assertthat::assert_that(is.numeric(mutations[,alt_reads_column]), 
    	all(mutations[,alt_reads_column] %% 1 == 0), all(mutations[,alt_reads_column] %% 1 == 0),
    	msg = "specified alt_reads_column should have all integer numbers")

    assertthat::assert_that(is.numeric(fdr_thr), length(fdr_thr) == 1,
        fdr_thr <= 1, fdr_thr > 0)

    if(!is.na(cfDNA_conc)){

    	assertthat::assert_that(is.numeric(cfDNA_conc), length(cfDNA_conc) == 1)
    
    }
    
    results <- list(VAF = NA, ctDNA_level = NA, clustering = NA, 
    	selected_cluster = NA, cluster_chr_pvalue = NA, mutations = NA)
    
    if(mean(mutations[,vaf_column]) < 0.005 || median(mutations[,vaf_column]) == 0){
    	
    	warning("Vaf is too low. This function is intended to assess ctDNA level in samples with noticeable levels of ctDNA", immediate.= T)
        
        return(results)
    
    }

    if(only_SNVs) {

    	mutations <- mutations %>%
    	   dplyr::filter(nchar(.data$REF) == 1 & nchar(.data$ALT) == 1)
    
        assertthat::assert_that(nrow(mutations) > 0, 
         	msg = "No entries in the mutations file after removing non-SNVs")
    
        results$mutations <- mutations
    
    }

    if(remove_traces) {
       
        if(!is.na(background_rate)){

       	    assertthat::assert_that(is.numeric(background_rate), length(background_rate) == 1, 
       	    	background_rate < 1, background_rate > 0)
       	   
       	    bg <- background_rate
       
        } else {
       
       	    assertthat::assert_that(!missing(bam), !missing(targets), !missing(reference))

            assertthat::assert_that(all(mutations$CHROM %in% GenomeInfoDb::seqnames(reference)), 
                all(targets$chr %in% GenomeInfoDb::seqnames(reference)), 
                msg = "Chromosomes in mutations and/or targets don't match the specified reference")

       	    assertthat::assert_that(is.character(bam), length(bam) == 1, file.exists(bam))
    
            assertthat::assert_that(is.character(tag), length(tag) == 1)
    
            assertthat::assert_that(all(get_bam_chr(bam) %in% GenomeInfoDb::seqnames(reference)), 
                msg = "Chromosomes in bam file don't match the specified reference")

       	    bg <- get_background_rate(bam = bam, targets = targets, 
       	        reference = reference, tag = tag, ...)$rate
       
        }

        p.binomial.test <- purrr::map2_dbl(mutations[,alt_reads_column], 
       	     mutations[,alt_reads_column] + mutations[,ref_reads_column],
       	     ~ binom.test(.x, .y, p = bg, alternative = "greater")$p.value)
    
        p.binomial.adjusted <- p.adjust(p.binomial.test, method = "fdr")

        mutations <- mutations[ p.binomial.adjusted < fdr_thr, ]
       
        message(sprintf("removed %s mutation traces", sum(p.binomial.adjusted >= fdr_thr)))

        assertthat::assert_that(nrow(mutations) > 0, 
          	msg = "Got out of mutations, try to relax the fdr_thr")

        if(nrow(mutations) < 5) {

            warning(sprintf("Number of mutations left (%s) is too low for reliable assessment", nrow(mutations)))
       
        }

        results$mutations <- mutations
    
    }

    if(use_clustering){

    	bic <- mclust::mclustBIC(mutations[,vaf_column, drop = F])

    	clustering <- mclust::Mclust(mutations[,vaf_column,drop = F], x = bic)

        mutations$mclust <- clustering$classification

        chr_per_cluster <- data.frame(CHROM = mutations$CHROM, mclust = mutations$mclust)
        
        ncluster <- length(clustering$parameters$mean)

        if(!is.na(override_selected_cluster) && override_selected_cluster %in% c(1:ncluster)) {

            selected_cluster <- override_selected_cluster
        
        } else {
          
            cluster_chisq_p <- purrr::map_dbl(c(1:ncluster), 
        	     ~ chisq.test(rbind(table(chr_per_cluster$CHROM), 
        	         table(chr_per_cluster[chr_per_cluster$mclust == .x, "CHROM"])))$p.value)
        
            selected_cluster <- which.max(cluster_chisq_p)
    
            results$cluster_chr_pvalue <- cluster_chisq_p
        
        }
        
        VAF <- as.numeric(clustering$parameters$mean[selected_cluster])

        results$clustering <- clustering
        results$selected_cluster <- selected_cluster
        results$mutations <- mutations
        
    } else {
       
       VAF <- mean(mutations[,vaf_column])

    }
    
    results$VAF <- VAF
    
    results$ctDNA_level <- VAF * cfDNA_conc * 3.3

    return(results)

}