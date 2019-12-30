#' Gets names of the reads showing reference and alternative allele of a list of mutations
#'
#' This function extracts the names of the reads in a bam file that support the variant and reference alleles of the input mutations
#' @param bam path to bam file
#' @param mutations A data frame containing the mutations. Must have the columns CHROM, POS, REF, ALT.
#' @param tag the RG tag if the bam has more than one sample
#' @param min_base_quality integer specifying the minimum base quality for reads to be included.
#' @param min_mapq integer specifying the minimum mapping quality for reads to be included
#' @return A list with length equal to the number of mutations. Each element is a character vector with the read names.
#' @seealso \code{\link{get_mutations_read_counts}} \code{\link{get_mutations_fragment_size}} \code{\link{test_ctDNA}}
#' @details Returns the IDs of the read that cover the input mutations (ref and alt alleles).
#' @export
#'
#' @examples
#' data("mutations", package = "ctDNAtools")
#' bamT1 <- system.file("extdata", "T1.bam", package = "ctDNAtools")
#' get_mutations_read_names(bam = bamT1, mutations = mutations[1:3, ])
get_mutations_read_names <- function(bam, mutations, min_base_quality = 20, tag = "", min_mapq = 30) {
  assertthat::assert_that(
    !missing(bam), is.character(bam),
    length(bam) == 1, file.exists(bam)
  )

  assertthat::assert_that(
    !missing(mutations),
    is.data.frame(mutations), assertthat::not_empty(mutations),
    assertthat::has_name(mutations, c("CHROM", "POS", "REF", "ALT"))
  )

  assertthat::assert_that(all(nchar(mutations$REF) == 1),
    all(nchar(mutations$ALT) == 1),
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

  read_names <- purrr::pmap(
    list(mutations$CHROM, mutations$POS, mutations$REF, mutations$ALT),
    function(chr, pos, ref, alt) {
      counts <- get_mutation_read_names(
        bam = bam, tag = tag, chr = chr, pos = pos,
        ref = ref, alt = alt, min_base_quality = min_base_quality, min_mapq = min_mapq
      )
    }
  )

  names(read_names) <- paste0(mutations$CHROM, ":", mutations$POS, "_", mutations$REF, "_", mutations$ALT)

  return(read_names)
}
