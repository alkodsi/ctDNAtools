#' Helper function to read a vcf into the required format of mutations data frame
#'
#' Uses VariantAnnotation::readVcfAsVRanges to read the vcf file, which return variants
#' in a format that each row is one variant. If the vcf has multiple samples, the samples
#' will be appended by rows. Provide a sample_name to return only the variants belonging
#' to the sample of interest. Once you use this function, make sure that all the variants
#' are relevant. The function will only return SNVs.

#' @param vcf the path to vcf file
#' @param sample_name a character(1) when provided, return only variants from this sample
#' @param ... other options passed to VariantAnnotation::readVcfAsVRanges
#'
#' @export
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#'
#' @examples
#' \donttest{
#' vcf <- system.file("extdata", "chr22.vcf.gz", package="VariantAnnotation")
#' vcf_to_mutations_df(vcf, sample_name = "HG00096")
#' }
#'
vcf_to_mutations_df <- function(vcf, sample_name = NULL, ...) {
  assertthat::assert_that(is.character(vcf), length(vcf) == 1, file.exists(vcf))

  ellipsis::check_dots_used()

  mutations <- as.data.frame(VariantAnnotation::readVcfAsVRanges(vcf, ...)) %>%
    dplyr::filter(nchar(.data$ref) == 1 & nchar(.data$alt) == 1) %>%
    dplyr::select(
      CHROM = .data$seqnames, POS = .data$start,
      REF = .data$ref, ALT = .data$alt, dplyr::everything()
    ) %>%
    dplyr::mutate(CHROM = as.character(.data$CHROM))

  if (!is.null(sample_name)) {
    assertthat::assert_that(is.character(sample_name), length(sample_name) == 1)

    assertthat::assert_that(sample_name %in% mutations$sampleNames,
      msg = "specified sample_name is not within the vcf"
    )

    mutations <- mutations %>%
      dplyr::filter(.data$sampleNames == sample_name)
  }

  return(mutations)
}
