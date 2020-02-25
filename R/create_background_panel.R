
#' Creates a background panel from a list of bam files
#'
#' This function scans the targets regions in the list of bam files, and reports the number of reference, non-reference reads for each loci
#' in addition to the non-reference (VAF) allele frequency. Loci with VAF higher than vaf_threshold are masked with NA.
#' @param bam_list A character vector of paths to bam files.
#' @param targets The targets data frame must have the columns chr, start and end.
#' @param reference The reference genome as BSgenome object.
#' @param vaf_threshold Loci with the fraction of non-reference reads above this value are masked with NA.
#' @param bam_list_tags RG tags for the list of bam files. By default, the whole bam file will be used.
#' @param min_base_quality The minimum base quality to count a read for a loci.
#' @param max_depth Maximum depth for the pileup
#' @param min_mapq The minimum mapping quality to count a read for a loci
#' @param substitution_specific logical, whether to have the loci by substitutions.
#' @return A named list having depth, alt and vaf data frames. Each has the same order of loci in rows and the input samples in columns.
#' @seealso  \code{\link{create_black_list}}  \code{\link{test_ctDNA}}
#' @details Extracts the depth, variant allele counts and variant allele frequency (VAF) for each genomic position in the input targets
#' across a panel of bam files (e.g. from healthy samples to represent only technical noise). The extracted information can be fed to
#' \code{\link{create_black_list}} in order to extract a black listed loci according to defined criteria
#'
#' The function support two modes, either loci-specific regardless of the basepair substitution, or substitution-specific where each
#' substitution class (e.g. C>T, C>G) are quantified separately. This behavior is controlled by the substitution_specific parameter.
#'
#' VAF above vaf_threshold parameters are masked with NA, to exclude real SNPs/mutations.
#'
#' Since this function can take a long time when the bam_list comprises a large number of bam files, the function supports multi-threading
#' using the furrr and future R packages. All you need to do is call 'plan(multiprocess)' or other multi-threading strategies before
#' calling this function.
#' @export
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
#'   min_samples_one_read = 2, min_samples_two_reads = 1
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
#'   min_samples_one_read = 2, min_samples_two_reads = 1
#' )
#'
#' ## Multi-threading (commented)
#' ## library(furrr)
#' ## plan(multiprocess)
#' ## plan(multiprocess, workers = 3)
#' bg_panel <- create_background_panel(
#'   bam_list = c(bamN1, bamN2, bamN3),
#'   targets = targets, reference = BSgenome.Hsapiens.UCSC.hg19,
#'   substitution_specific = TRUE
#' )
#' }
#'
create_background_panel <- function(bam_list, targets, reference, vaf_threshold = 0.05, bam_list_tags = rep("", length(bam_list)),
                                    min_base_quality = 10, max_depth = 1e+05, min_mapq = 20, substitution_specific = TRUE) {
  assertthat::assert_that(class(reference) == "BSgenome")
  assertthat::assert_that(
    is.data.frame(targets), assertthat::not_empty(targets),
    assertthat::has_name(targets, c("chr", "start", "end"))
  )

  assertthat::assert_that(is.character(bam_list), all(file.exists(bam_list)))

  assertthat::assert_that(length(bam_list) > 1,
    msg = "At least two bam files should be provided"
  )

  assertthat::assert_that(is.logical(substitution_specific) & length(substitution_specific) == 1)

  bam_list_chr <- purrr::map(bam_list, get_bam_chr)

  assertthat::assert_that(all(purrr::map_lgl(bam_list_chr, ~ all(.x %in% GenomeInfoDb::seqnames(reference)))),
    msg = "The chromosomes in at least one of the specified bams in bam_list don't match the reference"
  )

  assertthat::assert_that(
    is.character(bam_list_tags),
    length(bam_list_tags) == length(bam_list)
  )

  if (any(bam_list_tags != "")) {
    tag_verification <- purrr::map2_lgl(
      bam_list[bam_list_tags != ""],
      bam_list_tags[bam_list_tags != ""], verify_tag
    )

    assertthat::assert_that(all(tag_verification),
      msg = paste("specified tag for", paste(bam_list[bam_list_tags != ""][!tag_verification], collapse = " , "), "is not correct")
    )
  }

  background_panel <- furrr::future_map2(bam_list, bam_list_tags,
    ~ create_background_panel_instance(
      bam = .x, tag = .y, targets = targets, reference = reference,
      vaf_threshold = vaf_threshold, min_base_quality = min_base_quality, min_mapq = min_mapq,
      max_depth = max_depth, substitution_specific = substitution_specific
    )
  )

  sm <- make.unique(purrr::map_chr(bam_list, get_bam_SM))

  depth <- purrr::reduce(purrr::map(background_panel, "depth"), dplyr::full_join, by = "Locus") %>%
    rlang::set_names(c("Locus", sm))

  alt <- purrr::reduce(purrr::map(background_panel, "alt"), dplyr::full_join, by = "Locus") %>%
    rlang::set_names(c("Locus", sm))

  vaf <- purrr::reduce(purrr::map(background_panel, "vaf"), dplyr::full_join, by = "Locus") %>%
    rlang::set_names(c("Locus", sm))

  return(list(depth = depth, alt = alt, vaf = vaf))
}
