# ` Helper function
#'
#' Compares simulated and observed variant allele counts
#' @param simulated the simulated variant allele counts
#' @param observed the observed variant allele counts
#' @param depths the total depths
#' @keywords internal
compare_simulated_observed <- function(simulated, observed, depths) {
  observed_vaf <- mean(observed / depths, na.rm = TRUE)
  observed_nonzero <- sum(observed > 0)

  simulated_vaf <- mean(simulated / depths, na.rm = TRUE)
  simulated_nonzero <- sum(simulated > 0)

  out <- ifelse((simulated_vaf >= observed_vaf) &
    (simulated_nonzero >= observed_nonzero), 1, 0)
  return(out)
}


#' A function to sample from binomial distribution
#'
#' Samples from binomial distribution N variants with different depths and fixed mismatch rate for one simulation. Return 1 if the simulation exceeds observed alt allele reads or 0 otherwise
#' @param depths a vector with the depths of the variants
#' @param rate A list containing mismatch rates as produced by get_background_rate function
#' @param alt_reads the observed variant allele reads
#' @param seed the random seed
#' @return a scalar. Either 1 if the simulation exceeds observed variant alleles or 0 otherwise
#' @importFrom stats C rbinom runif
#' @seealso \code{\link{test_ctDNA}}
#' @keywords internal

simulator <- function(depths, rate, alt_reads, seed) {
  assertthat::assert_that(length(depths) == length(alt_reads))

  assertthat::assert_that(all(alt_reads <= depths))

  assertthat::assert_that(is.list(rate),
    assertthat::has_name(rate, c("rate", "CT", "CA", "CG", "TA", "TC", "TG")),
    msg = "rate must be a list as produced by get_background_rate"
  )

  n_variants <- length(depths)

  set.seed(seed)

  sim <- rbinom(n = n_variants, size = depths, prob = rate$rate / 3)
  out <- compare_simulated_observed(simulated = sim, observed = alt_reads, depths = depths)

  return(out)
}
