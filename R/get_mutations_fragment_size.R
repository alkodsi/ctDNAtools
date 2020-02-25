#' Gets reads fragment lengths for a list of mutations
#'
#' The function extracts the fragment lengths for the reads holding alternative allele for each mutation
#' in the mutations data frame.
#' @param bam path to bam file.
#' @param mutations Data frame with mutations. Must have the columns CHROM, POS, REF, ALT.
#' @param tag the RG tag if the bam has more than one sample.
#' @param min_base_quality minimum base quality when extracting reads covering mutations.
#' @param min_mapq minimum mapping quality when extracting reads covering mutations.
#' @param ... Other parameters passed to get_fragment_size.
#' @return A list with length equal to the number of mutations.
#' Each element contains a list with two elements ref and alt each having an integer vector of fragment lengths
#' @details Fragment length will extracted from the bam file according to the parameters passed to \code{\link{get_fragment_size}},
#' and the fragment size of the reads that map to the ref and alt alleles of each mutation in the input will be returned.
#' @seealso \code{\link{get_fragment_size}}
#' @export
#' @examples
#' \donttest{
#' data("mutations", package = "ctDNAtools")
#' bamT1 <- system.file("extdata", "T1.bam", package = "ctDNAtools")
#'
#' mfs <- get_mutations_fragment_size(bam = bamT1, mutations = mutations[1:2, ])
#' }
get_mutations_fragment_size <- function(bam, mutations, tag = "", min_base_quality = 20,
                                        min_mapq = 30, ...) {
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

  ellipsis::check_dots_used()

  mutations_targets <- data.frame(
    chr = mutations$CHROM, start = mutations$POS,
    end = mutations$POS, stringsAsFactors = FALSE
  )

  frag_size <- get_fragment_size(bam = bam, tag = tag, targets = mutations_targets, ...)

  read_names <- get_mutations_read_names(
    bam = bam, mutations = mutations, tag = tag,
    min_base_quality = min_base_quality, min_mapq = min_mapq
  )

  mutation_frag_size <- purrr::map(
    read_names,
    ~ list(ref = frag_size[frag_size$ID %in% .x$ref, "size"], alt = frag_size[frag_size$ID %in% .x$alt, "size"])
  )

  return(mutation_frag_size)
}
