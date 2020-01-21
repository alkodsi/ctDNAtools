#' Helper function to bin a variable

#' @param x the variable to be binned
#' @param from the minimum range
#' @param to the maximum range
#' @param by the break length
#' @param normalized whether to normalize the counts to total counts.
#' @param custom_bins A numeric vector for custom breaks to bin the histogram of fragment length. Over-rides bin_size.

#' @return A numeric vector having counts within bins normalized by the sum of the variable
#' @importFrom graphics hist
#' @keywords internal

get_hist_bins <- function(x, from, to, by, normalized, custom_bins = NULL) {
  assertthat::assert_that(is.numeric(x))

  if (!is.null(custom_bins)) {
    breaks <- unique(sort(c(from, custom_bins, to)))
  } else {
    breaks <- seq(from = from, to = to, by = by)
  }

  if (breaks[length(breaks)] != to) {
    breaks <- c(breaks, to)
  }

  h <- hist(x, breaks = breaks, plot = FALSE)

  if (normalized) {
    h$counts <- h$counts / sum(h$counts)
  }

  return(list(counts = h$counts, breaks = paste(h$breaks[-length(h$breaks)], h$breaks[-1], sep = "_")))
}
