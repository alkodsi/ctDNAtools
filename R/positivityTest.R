#' Compute a p-value of ctDNA positivity
#'
#' A function to determine ctDNA positivity with N repeated permutations. Return a p-value based on permutation test. Calls simulator function  nPermutation times.
#' @param depths a vector with the depths of the variants
#' @param altReads the observed variant allele reads
#' @param rate the mismatch rate
#' @param seed the random seed
#' @param nPermutation the number of permutations to run.
#' @return a scalar, permutation p-value.
#' @export

positivityTest <- function(depths, altReads, rate, seed, nPermutation) {
    set.seed(seed)
    
    assertthat::assert_that(assertthat::noNA(depths), assertthat::not_empty(depths))
    assertthat::assert_that(assertthat::noNA(altReads), assertthat::not_empty(altReads))
    
    seeds <- round(runif(n = nPermutation, min = 0, max = 1e+08))
    
    pvalue <- sum(purrr::map_dbl(seeds, ~ simulator(length(depths), depths, rate, 
        altReads, seed = .x)))/nPermutation
    
    return(pvalue)
}
