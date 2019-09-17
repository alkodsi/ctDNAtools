#' Compute a p-value of ctDNA positivity
#'
#' A function to determine ctDNA positivity with N repeated permutations. Return a p-value based on permutation test. Calls simulator function  n_permutation times.
#' @param depths a vector with the depths of the variants
#' @param altReads the observed variant allele reads
#' @param rate  A list containing mismatch rates as produced by get_background_rate function
#' @param substitutions character vector containing the substitutions e.g. CT, TG etc.
#' @param seed the random seed
#' @param n_permutation the number of permutations to run.
#' @return a scalar, permutation p-value.
#' @export

positivity_test <- function(depths, altReads, rate, substitutions = NULL, seed, n_permutation = 1e+04) {
    set.seed(seed)
    
    assertthat::assert_that(assertthat::noNA(depths), assertthat::not_empty(depths))
    assertthat::assert_that(assertthat::noNA(altReads), assertthat::not_empty(altReads))
    
    seeds <- round(runif(n = n_permutation, min = 0, max = 1e+08))
    
    pvalue <- sum(purrr::map_dbl(seeds, ~ simulator(depths = depths, rate = rate, 
        altReads = altReads, substitutions = substitutions, seed = .x)))/n_permutation
    
    return(pvalue)
}
