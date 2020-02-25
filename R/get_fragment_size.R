#' Gets fragment lengths from a bam file
#'
#' A function to extract fragment lengths from a bam file. Optionally, given a mutation data frame, it can categorize read lengths
#' in mutated vs non-mutated reads.
#' @param bam path to bam file.
#' @param mutations An optional data frame with mutations. Must have the columns CHROM, POS, REF, ALT.
#' @param tag the RG tag if the bam has more than one sample.
#' @param targets a data frame with the target regions to restrict the reads in the bam. Must have three columns: chr, start and end
#' @param isProperPair a logical whether to return only proper pairs (true), only improper pairs (false), or it does not matter (NA).
#' @param mapqFilter mapping quality threshold for considering reads.
#' @param min_size Integer with the lowest fragment length.
#' @param max_size Integer with the highest fragment length.
#' @param ignore_trimmed logical, whether to remove reads that have been hard trimmed.
#' @param different_strands logical, whether to keep only reads whose mates map to different strand.
#' @param simple_cigar logical, whether to include only reads with simple cigar.
#' @return A data frame with the columns:
#'  \itemize{
#'
#'    \item Sample: The SM tag in bam or file name
#'
#'    \item ID: the read ID
#'
#'    \item chr: chromosome
#'
#'    \item start: the left most end of either the read or mate
#'
#'    \item end: the right most end of either the read or mate.
#'
#'    \item size:  the fragment size
#'
#'    \item category (only if mutations is provided): either ref, alt, or other
#'    }
#'
#' @export
#' @details Extracts the fragment size of reads in the input bam that satisfy the following conditions:
#'  \itemize{
#'   \item Paired, and optionally properly paired depending on the isProperPair parameter.
#'
#'   \item Both the reads and mate are mapped.
#'
#'   \item Not secondary or supplementary alignment
#'
#'   \item Not duplicate
#'
#'   \item Passing quality control
#'
#'   \item Read and mate on different strands (optional depending on the different_strands parameter)
#'
#'   \item Not trimmed (optional depending on the ignore_trimmed parameter), i.e. will keep only reads with the max length.
#'
#'   \item Having a simple cigar (optional depending on the simple_cigar parameter).
#'
#'   \item Satisfy the mapping quality threshold specified in the mapqFilter parameter.
#'
#'   \item Reads and mates on the same chromosome when min_size > 0.
#' }
#'
#' When the input mutations is given, the output will label the reads that support the variant alleles of the mutation in
#' a separate column.

#' @seealso \code{\link{summarize_fragment_size}} \code{\link{bin_fragment_size}} \code{\link{analyze_fragmentation}} \code{\link{get_mutations_fragment_size}}
#' @importFrom rlang .data
#' @importFrom magrittr %>%
#' @examples
#' \donttest{
#' data("mutations", package = "ctDNAtools")
#' data("targets", package = "ctDNAtools")
#' bamT1 <- system.file("extdata", "T1.bam", package = "ctDNAtools")
#'
#' ## basic usage
#' fs <- get_fragment_size(bam = bamT1)
#'
#' ## More options
#' fs1 <- get_fragment_size(
#'   bam = bamT1, isProperPair = TRUE, min_size = 70,
#'   max_size = 200, ignore_trimmed = FALSE, different_strands = FALSE,
#'   simple_cigar = TRUE
#' )
#'
#' ## with mutations input
#' fs2 <- get_fragment_size(bam = bamT1, mutations = mutations)
#'
#' ## using targets
#' fs3 <- get_fragment_size(bam = bamT1, targets = targets)
#' }
#'
get_fragment_size <- function(bam, mutations = NULL, targets = NULL, tag = "", isProperPair = NA, mapqFilter = 30,
                              min_size = 1, max_size = 400, ignore_trimmed = TRUE, different_strands = TRUE, simple_cigar = FALSE) {
  assertthat::assert_that(!missing(bam), is.character(bam), length(bam) == 1, file.exists(bam))

  assertthat::assert_that(is.character(tag), length(tag) == 1)

  assertthat::assert_that(
    is.logical(isProperPair), is.logical(ignore_trimmed),
    is.logical(different_strands), is.logical(simple_cigar),
    is.numeric(min_size), is.numeric(max_size), max_size > min_size,
    length(isProperPair) == 1, length(ignore_trimmed) == 1, length(different_strands) == 1,
    length(simple_cigar) == 1, length(min_size) == 1, length(max_size) == 1
  )

  if (!is.null(mutations)) {
    assertthat::assert_that(
      is.data.frame(mutations), assertthat::not_empty(mutations),

      assertthat::has_name(mutations, c("CHROM", "POS", "REF", "ALT"))
    )

    assertthat::assert_that(all(nchar(mutations$REF) == 1),
      all(nchar(mutations$ALT) == 1),
      msg = "Only SNVs are supported"
    )

    assertthat::assert_that(is.character(mutations$REF), is.character(mutations$ALT),
      all(mutations$REF %in% c("A", "C", "T", "G")), all(mutations$ALT %in% c("A", "C", "T", "G")),
      msg = "REF and ALT in mutations should be characters having basepairs"
    )

    assertthat::assert_that(
      !any(duplicated(mutations[,c("CHROM", "POS", "REF", "ALT")])),
      msg = "mutations input has duplicates")

    assertthat::assert_that(is.numeric(mutations$POS), all(mutations$POS > 0))

    assertthat::assert_that(all(mutations$CHROM %in% get_bam_chr(bam)))
  }

  flag <- Rsamtools::scanBamFlag(
    isPaired = TRUE, isProperPair = isProperPair, isUnmappedQuery = FALSE,
    hasUnmappedMate = FALSE, isSecondaryAlignment = FALSE,
    isNotPassingQualityControls = FALSE, isDuplicate = FALSE, isSupplementaryAlignment = FALSE
  )

  if (tag != "") {
    assertthat::assert_that(verify_tag(bam = bam, tag = tag),
      msg = "Specified tag not found"
    )

    sbp <- Rsamtools::ScanBamParam(
      flag = flag, mapqFilter = mapqFilter, simpleCigar = simple_cigar,
      tagFilter = list(RG = tag), what = c("qname", "flag", "qwidth", "isize")
    )
  } else {
    sbp <- Rsamtools::ScanBamParam(
      flag = flag, mapqFilter = mapqFilter, simpleCigar = simple_cigar,
      what = c("qname", "flag", "qwidth", "isize")
    )
  }

  if (!is.null(targets)) {
    assertthat::assert_that(
      assertthat::has_name(targets, c("chr", "start", "end")),
      is.data.frame(targets), assertthat::not_empty(targets)
    )

    assertthat::assert_that(
      is.numeric(targets$start), assertthat::noNA(targets$start),
      all(targets$start > 0), assertthat::noNA(targets$end),
      is.numeric(targets$end), all(targets$end > 0)
    )

    gr <- GenomicRanges::reduce(GenomicRanges::GRanges(targets$chr, IRanges::IRanges(targets$start, targets$end)))

    Rsamtools::bamWhich(sbp) <- gr
  }

  sm <- get_bam_SM(bam = bam, tag = tag)

  reads <- GenomicAlignments::readGAlignmentPairs(file = bam, param = sbp)

  reads <- reads[abs(reads@first@elementMetadata$isize) > min_size &
    abs(reads@first@elementMetadata$isize) < max_size]

  if (different_strands) {
    reads <- reads[xor(BiocGenerics::strand(reads@first) == "+",
      BiocGenerics::strand(reads@last) == "+")]
  }

  if (ignore_trimmed) {
    max_qwidth <- max(reads@first@elementMetadata$qwidth)

    reads <- reads[reads@first@elementMetadata$qwidth == max_qwidth &
      reads@last@elementMetadata$qwidth == max_qwidth]
  }

  fragment_lengths <- as.data.frame(reads) %>%
    dplyr::mutate(
      Sample = sm,
      ID = paste(sm, .data$qname.first, sep = "_"),
      chr = .data$seqnames.first,
      start = pmin(.data$start.first, .data$end.first, .data$start.last, .data$end.last),
      end = pmax(.data$start.first, .data$end.first, .data$start.last, .data$end.last),
      size = abs(.data$isize.first),
      stringsAsFactors = FALSE
    ) %>%
    dplyr::select(
      .data$Sample, .data$ID, .data$chr,
      .data$start, .data$end, .data$size
    ) %>%
    dplyr::filter(!duplicated(.data$ID)) # using targets may lead to dups

  if (!is.null(mutations)) {
    read_names <- get_mutations_read_names(bam = bam, tag = tag, mutations = mutations)

    reads_alt <- unique(unlist(purrr::map(read_names, "alt")))

    reads_ref <- unique(unlist(purrr::map(read_names, "ref")))

    fragment_lengths <- dplyr::mutate(fragment_lengths,
      category = ifelse(.data$ID %in% reads_alt, "alt",
        ifelse(.data$ID %in% reads_ref, "ref", "other")))
  }

  return(fragment_lengths)
}
