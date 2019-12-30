#' Collapses mutations in phase into one event
#'
#' Given a mutations data frame and a bam file, this function collapses mutations in phase identified by the ID_column into one event.
#' While doing that, it ignores the reads that support both the reference and alternative alleles for different mutations in phase.

#' @param mutations A data frame with the reporter mutations. Should have the columns CHROM, POS, REF, ALT.
#' @param bam path to bam file
#' @param tag the RG tag if the bam has more than one sample
#' @param ID_column The name of the column in mutations data.frame that has the IDs for mutations in phase.
#' NA values will be filled automatically by unique mutation identifiers.
#' @param min_base_quality minimum base quality for a read to be counted
#' @param min_mapq integer specifying the minimum mapping quality for reads to be included.
#' @return A list with the following slots:
#' \describe{
#'  \item{out:}{ A data frame that has the columns:
#'       \itemize{
#'         \item Phasing_id: the ID of the mutations/event.
#'         \item ref: number of reference reads.
#'         \item alt: number of alternative reads.
#'         \item n_reads_multi_mutation: Number of reads that span more than one mutation in phase.
#'         \item all_reads: total number of reads.
#'         \item multi_support: number of reads that support the alt allele of multiple mutations in phase.
#'       }}
#'  \item{purification_prob:}{Probability of purification: sum(n_reads_multi_mutation)/sum(all_reads)}
#'  \item{multi_support:}{Number of multi-support reads in all mutations/events}
#'  \item{informative_reads:}{Number of unique reads covering the mutations/events}
#' }
#' @export
#' @seealso \code{\link{test_ctDNA}} \code{\link{get_mutations_read_names}}
#' @details Mutations in phase are those that are supported by the same reads (same allele). The function doesn't identify mutations in phase, but rather use
#' an ID column in the input whose name is specified by ID_column to tell which mutations are in phase.
#'
#' Since two or more mutations can be supported by the same evidence,
#' this function merges these mutations into one event. The function will also remove the mismatches that are not exhibited in all the covered phased mutations (since
#' this function is developed for the intent of minimal residual disease testing).
#'
#' The output will include the merged mutations, the probability of purification, which is defined as the number of
#' reads covering at least two mutations in phase divided by the number of informative reads.
#' Informative reads count is the total number of unique reads mapping to the mutations input
#' (including both mutations in phase and other mutations).

#' @importFrom rlang .data
#' @importFrom magrittr %>%
#' @examples
#' data("mutations", package = "ctDNAtools")
#' bamT1 <- system.file("extdata", "T1.bam", package = "ctDNAtools")
#' merge_mutations_in_phase(mutations = mutations[5:10, ], bam = bamT1, ID_column = "PHASING")
merge_mutations_in_phase <- function(mutations, bam, tag = "", ID_column = "phasingID", min_base_quality = 20, min_mapq = 30) {
  assertthat::assert_that(
    is.data.frame(mutations), assertthat::not_empty(mutations),
    assertthat::has_name(mutations, c("CHROM", "POS", "REF", "ALT", ID_column))
  )

  assertthat::assert_that(
    !any(duplicated(mutations[,c("CHROM", "POS", "REF", "ALT")])),
    msg = "mutations input has duplicates")
  
  mutation_id <- paste0(mutations$CHROM, ":", mutations$POS, "_", mutations$REF, "_", mutations$ALT)
  mutations[, ID_column] <- as.character(ifelse(is.na(mutations[, ID_column]), mutation_id, mutations[, ID_column]))

  IDs <- data.frame(mutations_id = mutation_id, phasing_id = mutations[, ID_column], stringsAsFactors = FALSE)
  IDs_list <- purrr::map(unique(IDs$phasing_id), ~ IDs[IDs$phasing_id == .x, "mutations_id"])

  read_names <- get_mutations_read_names(
    mutations = mutations, bam = bam,
    tag = tag, min_base_quality = min_base_quality, min_mapq = min_mapq
  )

  out <- purrr::map(IDs_list, function(.x) {
    ## extract mutations belonging to the phase ID
    read_names_to_merge <- read_names[.x]

    ## get the ref and alt reads
    alt <- unlist(purrr::map(read_names_to_merge, "alt"))
    ref <- unlist(purrr::map(read_names_to_merge, "ref"))

    ## number of reads that map to more than one mutation
    all_reads <- c(ref, alt)
    n_reads_multi_mutation <- length(unique(all_reads[duplicated(all_reads)]))

    ## purify mismatches mapping only to one of the mutations in phase
    alt <- alt[!alt %in% ref]

    ## count how many reads support more than one mutation
    multi_support <- length(unique(alt[duplicated(alt)]))

    ## put the counts together
    df <- data.frame(
      ref = length(unique(ref)),
      alt = length(unique(alt)),
      n_reads_multi_mutation = n_reads_multi_mutation,
      all_reads = length(unique(all_reads)),
      multi_support = multi_support
    )
    ref <- unique(ref)
    alt <- unique(alt)
    list(df = df, ref = ref, alt = alt)
  })

  informative_reads <- length(unique(c(
    unlist(purrr::map(read_names, "ref")),
    unlist(purrr::map(read_names, "alt"))
  )))

  df <- purrr::map_dfr(out, "df")

  purification_prob <- sum(df$n_reads_multi_mutation) / sum(df$all_reads)

  out <- data.frame(
    Phasing_id = unique(IDs$phasing_id),
    df, stringsAsFactors = FALSE
  )

  return(list(
    out = out, purification_prob = purification_prob, multi_support = sum(df$multi_support),
    informative_reads = informative_reads
  ))
}
