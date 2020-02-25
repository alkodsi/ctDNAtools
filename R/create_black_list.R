
#' Creates a black list of genomic loci based on a background panel created from a list of bam files (e.g. healthy samples)
#'
#' The function applies criteria on the background panel to extract the noisy genomic loci. Criteria include minimum number of samples having
#' at least one, at least two, or at least n (\code{n_reads} parameter) non-reference allele. Additionally the quantile of mean VAF above which the loci are considered noisy

#' @param background_panel A list produced by create_background panel function
#' @param mean_vaf_quantile The quantile of mean VAF above which the loci are considered noisy. Use NA to skip this criterion.
#' @param min_samples_one_read Loci that at least this number of samples exhibit at least one non-reference reads are considered noisy. Use NA to skip this criterion.
#' @param min_samples_two_reads Loci that at least this number of samples exhibit at least two non-reference reads are considered noisy. Use NA to skip this criterion.
#' @param min_samples_n_reads Loci that at least this number of samples exhibit at least n non-reference reads (\code{n_reads} parameter) are considered noisy. Use NA to skip this criterion.
#' @param n_reads the number of reads to use in the \code{min_samples_n_reads} parameter
#' @return a character vector of the loci in the black list
#' @seealso \code{\link{create_background_panel}}  \code{\link{test_ctDNA}}
#' @export
#' @importFrom stats quantile
#'
#' @examples
#' \donttest{
#' ## Load example data
#' data("targets", package = "ctDNAtools")
#' bamN1 <- system.file("extdata", "N1.bam", package = "ctDNAtools")
#' bamN2 <- system.file("extdata", "N2.bam", package = "ctDNAtools")
#' bamN3 <- system.file("extdata", "N3.bam", package = "ctDNAtools")
#'
#' ## Use human reference genome from BSgenome.Hsapiens.UCSC.hg19 library
#' suppressMessages(library(BSgenome.Hsapiens.UCSC.hg19))
#'
#' ## Use a black list based on loci
#' bg_panel <- create_background_panel(
#'   bam_list = c(bamN1, bamN2, bamN3),
#'   targets = targets, reference = BSgenome.Hsapiens.UCSC.hg19,
#'   substitution_specific = FALSE
#' )
#'
#' bl1 <- create_black_list(bg_panel,
#'   mean_vaf_quantile = 0.99,
#'   min_samples_one_read = 2, min_samples_two_reads = 1,
#'   min_samples_n_reads = 1, n_reads = 3
#' )
#'
#' ## Use a substitution-specific black list
#' bg_panel <- create_background_panel(
#'   bam_list = c(bamN1, bamN2, bamN3),
#'   targets = targets, reference = BSgenome.Hsapiens.UCSC.hg19,
#'   substitution_specific = TRUE
#' )
#'
#' bl2 <- create_black_list(bg_panel,
#'   mean_vaf_quantile = 0.99,
#'   min_samples_one_read = 2, min_samples_two_reads = 1,
#'   min_samples_n_read = NA
#' )
#' }
#'
create_black_list <- function(background_panel, mean_vaf_quantile = 0.95,
                              min_samples_one_read = max(2, ceiling(ncol(background_panel$vaf) * 0.75)),
                              min_samples_two_reads = max(2, ceiling(ncol(background_panel$vaf) * 0.2)),
                              min_samples_n_reads = NA,
                              n_reads = NA) {
  assertthat::assert_that(is.list(background_panel),
    assertthat::has_name(background_panel, c("depth", "alt", "vaf")),
    msg = "background_panel should be a list as produced by create_background_panel"
  )

  if(!is.null(mean_vaf_quantile) && !is.na(mean_vaf_quantile)) {
    assertthat::assert_that(is.numeric(mean_vaf_quantile),
      length(mean_vaf_quantile) == 1, mean_vaf_quantile < 1, 
      mean_vaf_quantile > 0
    )

    mean_vaf <- rowMeans(background_panel$vaf[, -1, drop = FALSE], na.rm = TRUE)

    quant <- quantile(mean_vaf, mean_vaf_quantile, na.rm = TRUE)

    if (quant > 0) {
      mean_vaf_idx <- which(mean_vaf >= quantile(mean_vaf, mean_vaf_quantile, na.rm = TRUE))
      message(sprintf("%s loci added satisfying Mean VAF condition", length(mean_vaf_idx)))
    } else {
      mean_vaf_idx <- c()
      message("The quantile is zero, skipping this criterion")
    }
  } else {
    mean_vaf_idx <- c()
  }

  if(!is.null(min_samples_one_read) && !is.na(min_samples_one_read)) {
    assertthat::assert_that(is.numeric(min_samples_one_read),
      length(min_samples_one_read) == 1
    )
    
    loci_one_read <- rowSums(background_panel$alt[, -1, drop = FALSE] >= 1, na.rm = TRUE)
    loci_one_read_idx <- which(loci_one_read >= min_samples_one_read)
    message(sprintf("%s loci added satisfying one read condition", length(loci_one_read_idx)))
  } else {
    loci_one_read_idx <- c()
  }

  if(!is.null(min_samples_two_reads) && !is.na(min_samples_two_reads)) {
    assertthat::assert_that(is.numeric(min_samples_two_reads),
      length(min_samples_two_reads) == 1
    )
    
    loci_two_reads <- rowSums(background_panel$alt[, -1, drop = FALSE] >= 2, na.rm = TRUE)
    loci_two_reads_idx <- which(loci_two_reads >= min_samples_two_reads)
    message(sprintf("%s loci added satisfying two reads condition", length(loci_two_reads_idx)))
  } else {
    loci_two_reads_idx <- c()
  }

  if(!is.null(min_samples_n_reads) && !is.na(min_samples_n_reads) && !is.null(n_reads) && !is.na(n_reads)) {
    assertthat::assert_that(is.numeric(min_samples_n_reads),
      length(min_samples_n_reads) == 1
    )

    assertthat::assert_that(is.numeric(n_reads),
      length(n_reads) == 1
    )
   
    loci_n_reads <- rowSums(background_panel$alt[, -1, drop = FALSE] >= n_reads, na.rm = TRUE)
    loci_n_reads_idx <- which(loci_n_reads >= min_samples_n_reads)
    message(sprintf("%s loci added satisfying n = %s reads condition", length(loci_n_reads_idx), n_reads))
  } else {
    loci_n_reads_idx <- c()
  }

  idx <- unique(c(mean_vaf_idx, loci_one_read_idx, loci_two_reads_idx, loci_n_reads_idx))
  message(sprintf("Black list has %s loci", length(idx)))

  black_list <- background_panel$depth$Locus[idx]

  return(black_list)
}
