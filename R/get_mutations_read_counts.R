#' Counts ref and alt reads for a set of mutations
#'
#' Counts ref and alt reads for a set of mutations in a bam file
#' @param mutations A data frame with the reporter mutations. Should have the columns CHROM, POS, REF, ALT.
#' @param bam path to bam file
#' @param tag the RG tag if the bam has more than one sample
#' @param min_base_quality minimum base quality for a read to be counted
#' @param max_depth maximum depth above which sampling will happen
#' @param min_mapq the minimum mapping quality for a read to be counted
#' @return a named list contains: ref, vector of read counts of the reference alleles, and
#'         alt, vector of read counts of the alternative allele
#' @seealso \code{\link{get_mutations_read_names}} \code{\link{test_ctDNA}} \code{\link{get_mutations_fragment_size}}
#' @details Quantifies the reference and variant alleles for the input mutations in the input bam file. Useful for forced calling mutations.
#' @export
#'
#' @examples
#' data("mutations", package = "ctDNAtools")
#' bamT1 <- system.file("extdata", "T1.bam", package = "ctDNAtools")
#' get_mutations_read_counts(mutations = mutations[1:3, ], bam = bamT1)
get_mutations_read_counts <- function(mutations, bam, tag = "", min_base_quality = 20,
                                      max_depth = 1e+05, min_mapq = 30) {
  assertthat::assert_that(
    !missing(bam), is.character(bam),
    length(bam) == 1, file.exists(bam)
  )


  assertthat::assert_that(
    is.data.frame(mutations), assertthat::not_empty(mutations),
    assertthat::has_name(mutations, c("CHROM", "POS", "REF", "ALT"))
  )

  assertthat::assert_that(all(nchar(mutations$REF) == 1), all(nchar(mutations$ALT) == 1),
    msg = "Only SNVs are supported"
  )

  assertthat::assert_that(is.character(mutations$REF), is.character(mutations$ALT),
    all(mutations$REF %in% c("A", "C", "T", "G")),
    all(mutations$ALT %in% c("A", "C", "T", "G")),
    msg = "REF and ALT in mutations should be characters having basepairs"
  )

  assertthat::assert_that(
    !any(duplicated(mutations[,c("CHROM", "POS", "REF", "ALT")])),
    msg = "mutations input has duplicates")
  
  assertthat::assert_that(is.numeric(mutations$POS), all(mutations$POS > 0))

  assertthat::assert_that(all(mutations$CHROM %in% get_bam_chr(bam)))

  readcounts <- purrr::pmap_dfr(list(
    mutations$CHROM, mutations$POS, mutations$REF,
    mutations$ALT
  ), function(chr, pos, ref, alt) {
    counts <- get_read_counts(
      chr = chr, pos = pos, bam = bam, tag = tag, min_base_quality = min_base_quality,
      min_mapq = min_mapq, max_depth = max_depth
    )
    data.frame(ref = counts[[ref]], alt = counts[[alt]])
  })

  return(list(ref = readcounts$ref, alt = readcounts$alt))
}
