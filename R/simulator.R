#` Helper function
#'
#' Compares simulated and observed
compare_simulated_observed <- function(simulated, observed){

	  comparison <- purrr::map_dbl(c(0:max(observed)),
       ~ as.numeric(sum(simulated >= .x) >= sum(observed >= .x)))
    
      out <- ifelse(sum(comparison) == length(comparison), 1, 0)
      return(out)

}


#' A function to sample from binomial distribution
#'
#' Samples from binomial distribution N variants with different depths and fixed mismatch rate for one simulation. Return 1 if the simulation exceeds observed alt allele reads or 0 otherwise
#' @param depths a vector with the depths of the variants
#' @param rate A list containing mismatch rates as produced by get_background_rate function
#' @param altReads the observed variant allele reads
#' @param substitutions character vector containing the substitutions e.g. CT, TG etc.
#' @param seed the random seed
#' @return a scalar. Either 1 if the simulation exceeds observed variant alleles or 0 otherwise

simulator <- function(depths, rate, altReads, substitutions = NULL, seed) {
    
    assertthat::assert_that(length(depths) == length(altReads))
    assertthat::assert_that(is.list(rate), 
    	           assertthat::has_name(rate, c("rate", "CT", "CA", "CG","TA","TC","TG")),
    	           msg = "rate must be a list as produced by get_background_rate")

    n_variants <- length(depths)

    set.seed(seed)
    
    if(is.null(substitutions)){
    
      sim <- rbinom(n = n_variants, size = depths, prob = rate$rate/3)
      out <- compare_simulated_observed(simulated = sim, observed = altReads)

    } else {
    
      assertthat::assert_that(length(depths) == length(substitutions))

      depths_list <- split(depths, substitutions)
      altReads_list <- split(altReads, substitutions)
    
      sim <- map(names(depths_list),
      	  ~ rbinom(n = length(depths_list[[.x]]), size = depths_list[[.x]], prob = rate[[.x]]))

      comparison_by_sub <- map2_dbl(sim, altReads_list,
      	  ~ compare_simulated_observed(simulated = .x, observed = .y))
      out <- ifelse(sum(comparison_by_sub) == length(comparison_by_sub), 1, 0)
    }
    
    return(out)
}
