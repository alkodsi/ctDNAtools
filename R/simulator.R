#' A function to sample from binomial distribution
#'
#' Samples from binomial distribution N variants with different depths and fixed mismatch rate for one permutation. Return 1 if the simulation exceeds observed alt allele reads or 0 otherwise
#' @param nVariants the number of variants to simulate.
#' @param depths a vector with the depths of the variants
#' @param rate the mismatch rate
#' @param altReads the observed variant allele reads
#' @param seed the random seed
#' @return a scalar. Either 1 if the simulation exceedds observed variant alleles or 0 otherwise

simulator <- function(nVariants, depths, rate, altReads, seed) {
    
    set.seed(seed)
    
    sim <- rbinom(n = nVariants, size = depths, prob = rate)
    
    comparison <- purrr::map_dbl(c(0:max(altReads)),
       ~as.numeric(sum(sim >= .x) >= sum(altReads >= .x)))
    
    out <- ifelse(sum(comparison) == length(comparison), 1, 0)
    
    return(out)
}
