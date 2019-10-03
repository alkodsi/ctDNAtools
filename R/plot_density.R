#' Plot density of variable
#' @param df data frame in long format
#' @param var the column of which to plot the histogram
#' @param color the column containing the categories
#' @param binwidth bin width of the histogram
#' @return a ggplot object
#' @export

plot_density <- function(df, var = "size", color = "Sample", binwidth = 2.5) {
 
    ggplot2::ggplot() + 
    ggplot2::geom_freqpoly(data = df, 
  	   ggplot2::aes_string(x = var, y = "..density..", color = color), 
  	   binwidth = binwidth) +
    ggplot2::theme_minimal()

}