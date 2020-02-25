#' Tests the ctDNA positivity of a sample
#'
#' Given a set of reporter mutation, this functions counts the reads matching the reporter mutations in the sample to be tested,
#' estimates the mismatch rate for the sample to be tested, and then runs a Monte Carlo simulation test to determine whether the tested sample is positive or negative.
#' @param mutations A data frame with the reporter mutations. Should have the columns CHROM, POS, REF, ALT.
#' @param bam path to bam file
#' @param targets a data frame with the target regions. Must have three columns: chr, start and end
#' @param reference the reference genome in BSgenome format
#' @param vaf_threshold When calculating the background rate, the bases with higher than this VAF threshold will be ignored (real mutations/SNPs). 
#' @param tag the RG tag if the bam has more than one sample
#' @param ID_column The name of the column that contains the ID of mutations in phase. All mutations in Phase should have the same ID in that column.
#' @param black_list a character vector of genomic loci to filter. The format is chr_pos if substitution_specific is false or
#' chr_pos_ref_alt if substitution_specific is true. If given, will override \code{bam_list}.
#' @param substitution_specific logical, whether to have the loci of black_list by substitutions.
#' @param min_base_quality minimum base quality for a read to be counted
#' @param max_depth maximum depth above which sampling will happen
#' @param min_mapq the minimum mapping quality for a read to be counted
#' @param bam_list A vector containing the paths to bam files used to filter mutations.
#' Mutations that have more than min_alt_reads in more than min_samples will be filtered. Using black_list is more recommended.
#' @param bam_list_tags the RG tags for the bams included in bams list.
#' @param min_alt_reads When bam_list is provided, this sets the minimum number of alternative allele reads for a sample to be counted.
#' @param min_samples When number of samples having more than \code{min_alt_reads} exceeds this number, the mutation will be filtered.
#' @param n_simulations the number of Monte Carlo simulations.
#' @param pvalue_threshold the p-value threshold used to decide positivity or negativity.
#' @param seed the random seed to make the Monte Carlo simulations reproducible.
#' @param informative_reads_threshold the number of informative reads (unique reads mapping to specified mutations) under which the test will be undetermined.
#' @return a data frame with the following columns:
#' \itemize{
#'   \item sample: The sample name taken from SM field in the bam file or file base name
#'   \item n_mutations: The number of mutations used in the test.
#'   \item n_nonzero_alt: Number of mutations that have at least one read supporting alternative allele.
#'   \item total_alt_reads: Total number of reads supporting alternative alleles of all mutations in input.
#'   \item mutations_filtered: The number of filtered mutations.
#'   \item background_rate: The background rate of the tested sample (after all adjustments)
#'   \item informative_reads: The number of unique reads covering the mutations used.
#'   \item multi_support_reads: The number of reads that support more than one mutations in phase. Non-zero values is a sign of
#'    positivity not used in the p-value calculation.
#'   \item pvalue: The empirical p-value from the Monte Carlo test.
#'   \item decision: The decision can be positive, negative or undetermined.
#' }
#'
#' @export
#' @seealso \code{\link{get_background_rate}} \code{\link{merge_mutations_in_phase}} \code{\link{create_black_list}} \code{\link{create_background_panel}}
#'    \code{\link{filter_mutations}}

#' @details This is the main function to test minimal residual disease by ctDNA positivity in a follow-up sample (e.g. after treatment).
#' The inputs include a bam file for the follow-up sample to be tested, a list of reporter mutations (detected for example before treatment in a ctDNA positive sample),
#' and an optional black_list (recommended to use) computed from a list of bam files of healthy-like samples or bam_list of the bam_files to use instead of black_list.
#'
#' The workflow includes the following steps:
#'
#' \describe{
#' \item{1.}{Filtering mutations (optional but recommended): The mutations in the input will be filtered removing the ones reported in the black list. If bam_list is provided,
#'    the mutations will be filtered according to the min_samples and min_alt_reads parameters. See \code{\link{filter_mutations}}.}
#'
#' \item{2.}{The background rate will be computed for the input bam. The black list will be plugged in to exclude the black-listed loci when computing the background rate.
#'    The black list can be substitution_specific or not, but in both cases, the background rate will be adjusted accordingly. See \code{\link{get_background_rate}}.}
#'
#' \item{3.}{Counting reference and alternative alleles of the reporter mutations in the bam file.}
#'
#' \item{4.}{Merging mutations in phase (optional but recommended): If the ID_column is specified in the mutations input, mutations with the same ID (in phase) will be merged.
#'    While doing so, it is expected that real traces of mutations will be exhibited jointly in the reads spanning phased mutations, whereas artifacts will not show in all the covered
#'    mutations in phase. Therefore, the mismatches that map only to a subset of the mutations in phase but not the others (which are covered and show reference allele) will be removed.
#'    This will lead to reduction of the observed mismatches, and therefore, the computed background rate will be adjusted accordingly:
#'    new rate = old rate * (1 - purification_probability). See \code{\link{merge_mutations_in_phase}}.}
#'
#' \item{5.}{Monte Carlo test: this approach was used by Newman et al., Nature Biotechnology 2016.
#'    \itemize{
#'    \item  Given \eqn{N} reporter mutations each with depth \ifelse{html}{\out{D<sub>i</sub>}}{\eqn{D_i}}, randomly sample variant allele reads \ifelse{html}{\out{X<sub>i</sub>}}{\eqn{X_i}}
#'    under the background rate \eqn{p} using binomial distribution \ifelse{html}{\out{X<sub>i</sub>}}{\eqn{X_i}} ~ Binom(\ifelse{html}{\out{D<sub>i</sub>}}{\eqn{D_i}}, \eqn{p}).
#'
#'    \item  Repeat \code{n_simuations} times.
#'
#'    \item  Count the number of simulations, where simulated data equal or exceed observed data in jointly two measurements: (1) the average VAF for the \eqn{N} mutations,
#'      and (2) the number of mutations with non-zero VAF.
#'
#'    \item  Compute an empirical p-value as (#successes + 1)/(#simulations + 1)
#'    }}
#'
#' \item{6.}{Make a decision: If number of informative reads is lower than \code{informative_reads_threshold} parameter,
#'  the decision will be undetermined. Otherwise, the \code{pvalue_threshold}
#'    parameters will be used to determine positivity or negativity.}
#' }
#'
#' @examples
#' ## Load example data
#' data("mutations", package = "ctDNAtools")
#' data("targets", package = "ctDNAtools")
#' bamT1 <- system.file("extdata", "T1.bam", package = "ctDNAtools")
#' bamT2 <- system.file("extdata", "T2.bam", package = "ctDNAtools")
#' bamN1 <- system.file("extdata", "N1.bam", package = "ctDNAtools")
#' bamN2 <- system.file("extdata", "N2.bam", package = "ctDNAtools")
#' bamN3 <- system.file("extdata", "N3.bam", package = "ctDNAtools")
#'
#' ## Use human reference genome from BSgenome.Hsapiens.UCSC.hg19 library
#' suppressMessages(library(BSgenome.Hsapiens.UCSC.hg19))
#' \donttest{
#' ## basic usage
#' test_ctDNA(
#'   mutations = mutations, bam = bamT1,
#'   targets = targets, reference = BSgenome.Hsapiens.UCSC.hg19,
#'   n_simulation = 100
#' )
#'
#' ## More options
#' test_ctDNA(
#'   mutations = mutations, bam = bamT1,
#'   targets = targets, reference = BSgenome.Hsapiens.UCSC.hg19,
#'   n_simulation = 100, min_base_quality = 20, min_mapq = 30,
#'   vaf_threshold = 0.05
#' )
#'
#' ## Use phasing information
#' test_ctDNA(
#'   mutations = mutations, bam = bamT2,
#'   targets = targets, reference = BSgenome.Hsapiens.UCSC.hg19,
#'   n_simulation = 100, ID_column = "PHASING"
#' )
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
#' test_ctDNA(
#'   mutations = mutations, bam = bamT1,
#'   targets = targets, reference = BSgenome.Hsapiens.UCSC.hg19,
#'   n_simulation = 100, ID_column = "PHASING", black_list = bl1,
#'   substitution_specific = FALSE
#' )
#' }
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
#' test_ctDNA(
#'   mutations = mutations, bam = bamT1,
#'   targets = targets, reference = BSgenome.Hsapiens.UCSC.hg19,
#'   n_simulation = 100, ID_column = "PHASING", black_list = bl2,
#'   substitution_specific = TRUE
#' )
test_ctDNA <- function(mutations, bam, targets, reference, tag = "", ID_column = NULL, black_list = NULL, substitution_specific = TRUE,
                       vaf_threshold = 0.1, min_base_quality = 30, max_depth = 1e+05, min_mapq = 40, bam_list = NULL,
                       bam_list_tags = rep("", length(bam_list)), min_alt_reads = 1, min_samples = ceiling(length(bam_list) / 10),
                       n_simulations = 10000, pvalue_threshold = 0.05, seed = 123,
                       informative_reads_threshold = 10000) {
  assertthat::assert_that(!missing(mutations), !missing(bam), !missing(targets),
    !missing(reference),
    msg = "mutations, bam, targets and reference are all required"
  )

  assertthat::assert_that(
    is.data.frame(mutations),
    assertthat::has_name(mutations, c("CHROM", "POS", "REF", "ALT"))
  )

  assertthat::assert_that(
    is.data.frame(targets), assertthat::not_empty(targets),
    assertthat::has_name(targets, c("chr", "start", "end"))
  )

  assertthat::assert_that(class(reference) == "BSgenome")

  assertthat::assert_that(all(nchar(mutations$REF) == 1), all(nchar(mutations$ALT) == 1),
    msg = "Only SNVs are supported"
  )

  assertthat::assert_that(is.character(mutations$REF), is.character(mutations$ALT),
    all(mutations$REF %in% c("A", "C", "T", "G")), all(mutations$ALT %in% c("A", "C", "T", "G")),
    msg = "REF and ALT in mutations should be characters having basepairs"
  )

  assertthat::assert_that(
    is.numeric(mutations$POS), assertthat::noNA(mutations$POS),
    all(mutations$POS > 0)
  )

  assertthat::assert_that(
    !any(duplicated(mutations[,c("CHROM", "POS", "REF", "ALT")])),
    msg = "mutations input has duplicates")
  
  assertthat::assert_that(
    is.numeric(targets$start), assertthat::noNA(targets$start),
    all(targets$start > 0), assertthat::noNA(targets$end),
    is.numeric(targets$end), all(targets$end > 0)
  )

  assertthat::assert_that(all(mutations$CHROM %in% GenomeInfoDb::seqnames(reference)),
    all(targets$chr %in% GenomeInfoDb::seqnames(reference)),
    msg = "Chromosomes in mutations and/or targets don't match the specified reference"
  )

  assertthat::assert_that(is.character(bam), length(bam) == 1, file.exists(bam))

  assertthat::assert_that(is.character(tag), length(tag) == 1)

  assertthat::assert_that(all(get_bam_chr(bam) %in% GenomeInfoDb::seqnames(reference)),
    msg = "Chromosomes in bam file don't match the specified reference"
  )

  assertthat::assert_that(
    is.numeric(vaf_threshold), vaf_threshold > 0,
    vaf_threshold <= 1, length(vaf_threshold) == 1
  )

  assertthat::assert_that(
    is.numeric(min_base_quality), length(min_base_quality) == 1,
    min_base_quality %% 1 == 0, min_base_quality > 0
  )

  assertthat::assert_that(
    is.numeric(max_depth), length(max_depth) == 1,
    max_depth %% 1 == 0, max_depth > 0
  )

  assertthat::assert_that(
    is.numeric(min_mapq), length(min_mapq) == 1,
    min_mapq %% 1 == 0, min_mapq > 0
  )

  assertthat::assert_that(
    is.numeric(n_simulations), length(n_simulations) == 1,
    n_simulations %% 1 == 0, n_simulations > 0
  )

  assertthat::assert_that(
    is.numeric(seed), length(seed) == 1, seed %% 1 == 0,
    seed >= 0
  )

  if (tag != "") {
    assertthat::assert_that(verify_tag(bam = bam, tag = tag), msg = "Specified tag not found")
  }

  if (!is.null(ID_column)) {
    assertthat::assert_that(assertthat::has_name(mutations, ID_column))

    assertthat::assert_that(is.character(mutations[, ID_column]))
  }

  if (!is.null(bam_list)) {
    assertthat::assert_that(
      is.numeric(min_alt_reads), length(min_alt_reads) == 1,
      min_alt_reads %% 1 == 0, min_alt_reads >= 0
    )

    assertthat::assert_that(
      is.numeric(min_samples), length(min_samples) == 1,
      min_samples %% 1 == 0, min_samples >= 0
    )

    assertthat::assert_that(is.character(bam_list), all(file.exists(bam_list)))

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
  }

  if (!is.null(black_list)) {
    assertthat::assert_that(is.character(black_list))

    if (substitution_specific) {
      assertthat::assert_that(all(purrr::map_dbl(strsplit(black_list, "_"), length) == 4),
        all(purrr::map_chr(strsplit(black_list, "_"), 3) %in% c("C", "A", "T", "G")),
        all(purrr::map_chr(strsplit(black_list, "_"), 4) %in% c("C", "A", "T", "G")),
        msg = "black_list must have characters in the format chr_pos_ref_alt when substitution_specific is true"
      )
    } else {
      assertthat::assert_that(all(purrr::map_dbl(strsplit(black_list, "_"), length) == 2),
        msg = "black_list should have characters in the format chr_pos"
      )
    }


    assertthat::assert_that(all(purrr::map_chr(strsplit(black_list, "_"), 1) %in% GenomeInfoDb::seqnames(reference)),
      msg = "Chromosomes of black_list are not in reference"
    )
  }

  sm <- get_bam_SM(bam = bam, tag = tag)

  message(paste("Analyzing Sample", sm, "..."))

  if (nrow(mutations) == 0) {
    warning(sprintf("In Sample %s, no mutations in input", sm), immediate. = TRUE)

    out <- data.frame(
      sample = sm,
      n_mutations = nrow(mutations),
      n_nonzero_alt = NA,
      total_alt_reads = NA,
      mutations_filtered = NA,
      background_rate = NA,
      informative_reads = NA,
      multi_support_reads = NA,
      pvalue = NA,
      decision = factor("undetermined", levels = c("positive", "negative", "undetermined")),
      stringsAsFactors = FALSE
    )

    return(out)
  }

  n_filtered <- 0

  if (!is.null(bam_list) || !is.null(black_list)) {
    original_n <- nrow(mutations)

    mutations <- filter_mutations(
      mutations = mutations, bams = bam_list, black_list = black_list, tags = bam_list_tags,
      min_alt_reads = min_alt_reads, min_samples = min_samples, min_base_quality = min_base_quality,
      max_depth = max_depth, min_mapq = min_mapq, substitution_specific = substitution_specific
    )

    n_filtered <- original_n - nrow(mutations)

    if (nrow(mutations) == 0) {
      warning(sprintf("In Sample %s, all present mutations were filtered", sm), immediate. = TRUE)

      out <- data.frame(
        sample = sm,
        n_mutations = nrow(mutations),
        n_nonzero_alt = NA,
        total_alt_reads = NA,
        mutations_filtered = n_filtered,
        background_rate = NA,
        informative_reads = NA,
        multi_support_reads = NA,
        pvalue = NA,
        decision = factor("undetermined", levels = c("positive", "negative", "undetermined")),
        stringsAsFactors = FALSE
      )

      return(out)
    }
  }

  message("Estimating background rate ...")

  bg <- get_background_rate(
    bam = bam, targets = targets, reference = reference,
    tag = tag, vaf_threshold = vaf_threshold, min_base_quality = min_base_quality,
    black_list = black_list, max_depth = max_depth, min_mapq = min_mapq,
    substitution_specific = substitution_specific
  )

  message("Getting ref and alt Counts ...")

  if (!is.null(ID_column)) {
    message("merging mutations in phase ...")

    ref_alt_reads <- merge_mutations_in_phase(
      mutations = mutations, bam = bam,
      tag = tag, min_base_quality = min_base_quality, min_mapq = min_mapq, ID_column = ID_column
    )

    prob_purification <- ref_alt_reads$purification_prob

    informative_reads <- ref_alt_reads$informative_reads

    multi_support_reads <- ref_alt_reads$multi_support

    ref_alt_reads <- ref_alt_reads$out

    bg$rate <- bg$rate * (1 - prob_purification)

  } else {

    ref_alt_reads <- get_mutations_read_counts(
      mutations = mutations, bam = bam, tag = tag,
      min_base_quality = min_base_quality, max_depth = max_depth, min_mapq = min_mapq
    )

    read_names <- get_mutations_read_names(
      mutations = mutations, bam = bam,
      tag = tag, min_base_quality = min_base_quality, min_mapq = min_mapq
    )

    informative_reads <- length(unique(c(
      unlist(purrr::map(read_names, "ref")),
      unlist(purrr::map(read_names, "alt"))
    )))

    multi_support_reads <- NA
  }

  alt_reads <- ref_alt_reads$alt

  ref_reads <- ref_alt_reads$ref

  message("Running Monte Carlo simulations")

  pos_test <- positivity_test(
    depths = ref_reads + alt_reads, alt_reads = alt_reads,
    rate = bg, seed = seed, n_simulations = n_simulations
  )

  message(paste("Pvalue = ", pos_test))

  decision <- ifelse(pos_test < pvalue_threshold, "positive", "negative")

  decision <- ifelse(informative_reads < informative_reads_threshold, "undetermined", decision)

  message(paste("Sample", sm, "is", decision))

  out <- data.frame(
    sample = sm,
    n_mutations = nrow(mutations),
    n_nonzero_alt = sum(alt_reads > 0),
    total_alt_reads = sum(alt_reads),
    mutations_filtered = n_filtered,
    background_rate = bg$rate,
    informative_reads = informative_reads,
    multi_support_reads = multi_support_reads,
    pvalue = pos_test,
    decision = factor(decision, levels = c("positive", "negative", "undetermined")),
    stringsAsFactors = FALSE
  )

  return(out)
}
