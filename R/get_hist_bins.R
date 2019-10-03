#' Helper function to bin a variable

#' @param x the variable to be binned
#' @param from the minimum range
#' @param to the maximum range
#' @param by the break length
#' @return A numeric vector having counts within bins normalized by the sum of the variable

get_hist_bins <- function(x, from, to, by, normalized){
	
	assertthat::assert_that(is.numeric(x))

    breaks <- seq(from = from, to = to, by = by)

    if( breaks[length(breaks)] != to){

    	breaks <- c(breaks, to)
    
    }

    h <- hist(x, breaks = breaks, plot = F)
    
    if(normalized){
      
      return(h$counts/sum(h$counts))
    
    } else {
      
      return(h$counts)
    
    }
}