#' Computes a p-value of ctDNA positivity with a Monte-Carlo based sampling test
#'
#' A function to determine ctDNA positivity with N repeated Monte Carlo simulations. Return a p-value based on Monte Carlo simulation test. Calls simulator function  n_simulations times.
#' @param depths a vector with the depths of the variants
#' @param alt_reads the observed variant allele reads
#' @param rate  A named list containing mismatch rates as produced by get_background_rate function
#' @param seed the random seed
#' @param n_simulations the number of simulations to run.
#' @return a scalar, simulation p-value.
#' @seealso \code{\link{test_ctDNA}}
#' @keywords internal

positivity_test <- function(depths, alt_reads, rate, seed = 123,
                            n_simulations = 10000) {
  set.seed(seed)

  assertthat::assert_that(
    assertthat::noNA(depths), assertthat::not_empty(depths),
    is.numeric(depths)
  )

  assertthat::assert_that(
    assertthat::noNA(alt_reads), assertthat::not_empty(alt_reads),
    is.numeric(alt_reads)
  )

  assertthat::assert_that(is.numeric(depths), all(depths %% 1 == 0))

  assertthat::assert_that(is.numeric(alt_reads), all(alt_reads %% 1 == 0))

  assertthat::assert_that(
    length(depths) == length(alt_reads),
    all(depths - alt_reads >= 0)
  )

  assertthat::assert_that(
    is.numeric(n_simulations), length(n_simulations) == 1,
    n_simulations %% 1 == 0, n_simulations > 0
  )

  assertthat::assert_that(
    is.numeric(seed), length(seed) == 1, seed %% 1 == 0,
    seed >= 0
  )

  assertthat::assert_that(
    is.list(rate), assertthat::has_name(rate, "rate"),
    is.numeric(rate$rate), length(rate$rate) == 1
  )

  seeds <- round(runif(n = n_simulations, min = 0, max = 1e+08))

  pvalue <- (sum(purrr::map_dbl(seeds, ~ simulator(
    depths = depths, rate = rate,
    alt_reads = alt_reads, seed = .x))) + 1) / (n_simulations + 1)

  return(pvalue)
}
