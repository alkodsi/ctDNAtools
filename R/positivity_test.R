#' Compute a p-value of ctDNA positivity
#'
#' A function to determine ctDNA positivity with N repeated Monte Carlo simulations. Return a p-value based on Monte Carlo simulation test. Calls simulator function  n_simulations times.
#' @param depths a vector with the depths of the variants
#' @param altReads the observed variant allele reads
#' @param rate  A list containing mismatch rates as produced by get_background_rate function
#' @param substitutions character vector containing the substitutions e.g. CT, TG etc.
#' @param seed the random seed
#' @param n_simulations the number of simulations to run.
#' @return a scalar, simulation p-value.
#' @export

positivity_test <- function(depths, altReads, rate, substitutions = NULL, seed, n_simulations = 1e+04) {
    set.seed(seed)
    
    assertthat::assert_that(assertthat::noNA(depths), assertthat::not_empty(depths))
    assertthat::assert_that(assertthat::noNA(altReads), assertthat::not_empty(altReads))
    
    seeds <- round(runif(n = n_simulations, min = 0, max = 1e+08))
    
    pvalue <- (sum(purrr::map_dbl(seeds, ~ simulator(depths = depths, rate = rate, 
        altReads = altReads, substitutions = substitutions, seed = .x)))+1) / (n_simulations + 1)
    
    return(pvalue)
}
