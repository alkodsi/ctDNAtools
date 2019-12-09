
#' A function to create a black list of genomic loci based on a background panel created from a list of bam files (e.g. healthy samples)
#'
#' The function applies criteria on the background panel to extract the noisy genomic loci. Criteria include minimum number of samples having
#' at least one and at least two non-reference allele. Additionally the quantile of mean VAF above which the loci are considered noisy

#' @param background_panel A list produced by create_background panel function
#' @param mean_vaf_quantile The quantile of mean VAF above which the loci are considered noisy
#' @param min_samples_one_read Loci that at least this number of samples exhibit at least one non-reference reads are considered noisy.
#' @param min_samples_two_reads Loci that at least this number of samples exhibit at least two non-reference reads are considered noisy.
#' @return a character vector of the loci in the black list

#' @export
#' @importFrom stats quantile

create_black_list <- function(background_panel, mean_vaf_quantile = 0.95, 
    min_samples_one_read = max(2, ceiling(ncol(background_panel$vaf) * 0.75)), 
    min_samples_two_reads = max(2, ceiling(ncol(background_panel$vaf) * 0.2))) {

    assertthat::assert_that(is.list(background_panel), 
        assertthat::has_name(background_panel, c("depth","alt","vaf")),
        msg = "background_panel should be a list as produced by create_background_panel")
    
    mean_vaf <- rowMeans(background_panel$vaf[, -1, drop = F], na.rm = T)
    mean_vaf_idx <- which(mean_vaf >= quantile(mean_vaf, mean_vaf_quantile, na.rm = T))
    message(sprintf("%s loci added satisfying Mean VAF condition", length(mean_vaf_idx)))

    samples_one_read <- rowSums(background_panel$alt[, -1, drop = F] >= 1, na.rm =T)
    samples_one_read_idx <- which(samples_one_read >= min_samples_one_read)
    message(sprintf("%s loci added satisfying one read condition", length(samples_one_read_idx)))

    samples_two_reads <- rowSums(background_panel$alt[, -1, drop = F] >=2, na.rm = T)
    samples_two_reads_idx <- which(samples_two_reads >= min_samples_two_reads)
    message(sprintf("%s loci added satisfying two reads condition", length(samples_two_reads_idx)))

    idx <- unique(c(mean_vaf_idx, samples_one_read_idx, samples_two_reads_idx))
    message(sprintf("Black list has %s loci", length(idx)))

    black_list <- background_panel$depth$Locus[idx]
    
    return(black_list)
}