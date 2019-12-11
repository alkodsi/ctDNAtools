#' Compute a p-value of ctDNA positivity
#'
#' A function to determine ctDNA positivity with N repeated Monte Carlo simulations. Return a p-value based on Monte Carlo simulation test. Calls simulator function  n_simulations times.
#' @param depths a vector with the depths of the variants
#' @param altReads the observed variant allele reads
#' @param rate  A named list containing mismatch rates as produced by get_background_rate function
#' @param substitutions character vector containing the substitutions e.g. CT, TG etc.
#' @param seed the random seed
#' @param n_simulations the number of simulations to run.
#' @return a scalar, simulation p-value.

positivity_test <- function(depths, altReads, rate, substitutions = NULL, seed = 123, 
    n_simulations = 10000) {
    
    set.seed(seed)
    
    assertthat::assert_that(assertthat::noNA(depths), assertthat::not_empty(depths), 
        is.numeric(depths))
    
    assertthat::assert_that(assertthat::noNA(altReads), assertthat::not_empty(altReads), 
        is.numeric(altReads))
    
    assertthat::assert_that(is.numeric(depths), all(depths%%1 == 0))
    
    assertthat::assert_that(is.numeric(altReads), all(altReads%%1 == 0))
    
    assertthat::assert_that(length(depths) == length(altReads),
        all(depths - altReads >= 0))
    
    assertthat::assert_that(is.numeric(n_simulations), length(n_simulations) == 1, 
        n_simulations%%1 == 0, n_simulations > 0)
    
    assertthat::assert_that(is.numeric(seed), length(seed) == 1, seed%%1 == 0, 
        seed >= 0)
    
    assertthat::assert_that(is.list(rate), assertthat::has_name(rate, "rate"),
        is.numeric(rate$rate), length(rate$rate) == 1)
    
    if (!is.null(substitutions)) {
        assertthat::assert_that(is.character(substitutions),
            length(substitutions) == length(depths), all(nchar(substitutions) == 2), 
            all(substitutions %in% c("CT", "CA", "CG", "TA", "TC", "TG", "AT", "AG", "AC", "GC", "GA", "GT")))
        
        substitutions <- dplyr::case_when(
            substitutions %in% c("CT", "GA") ~ "CT", 
            substitutions %in% c("CA", "GT") ~ "CA",
            substitutions %in% c("CG", "GC") ~ "CG",
            substitutions %in% c("TA", "AT") ~ "TA", 
            substitutions %in% c("TC", "AG") ~ "TC",
            substitutions %in% c("TG", "AC") ~ "TG")
        
        assertthat::assert_that(assertthat::has_name(rate, c("CT", "CG", "CA", "TC", "TA", "TG")))
    }
    
    seeds <- round(runif(n = n_simulations, min = 0, max = 1e+08))
    
    pvalue <- (sum(purrr::map_dbl(seeds, ~simulator(depths = depths, rate = rate, 
        altReads = altReads, substitutions = substitutions, seed = .x))) + 1)/(n_simulations + 1)
    
    return(pvalue)
}
