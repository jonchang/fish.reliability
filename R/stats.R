#' 2D median absolute deviation
#'
#' Compute the median absolute deviation for two dimensional data.
#' @param x a numeric vector.
#' @param y a numeric vector, of the same length as \code{x}.
#' @param constant scale factor. See 'Details' in \code{\link[stats]{mad}}.
#' @export
#' @importFrom stats median
mad_2d <- function(x, y, constant = 1.4826) {
  if (length(x) != length(y))
    stop("'x' and 'y' must have the same length")
  cx <- median(x)
  cy <- median(y)

  constant * median(sqrt((x - cx)^2 + (y-cy)^2))
}

#' Median distance of two groups
#'
#' @param x a numeric vector
#' @param y a numeric vector
#' @param group a vector that specifies how \code{x} and \code{y} should be grouped
#' @export
#' @importFrom stats median dist
dist_2d <- function(x, y, group) {
  if (length(x) != length(y))
    stop("'x' and 'y' must have the same length")
  x <- as.numeric(x)
  y <- as.numeric(y)
  df <- data_frame(x, y, group) %>% group_by(group) %>% summarise(x=median(x), y=median(y)) %>% select(x,y)
  dist(df)[1]
}

#' Detects an upper (more positive) outlier
#'
#' @param x a numeric vector
#' @param y the cutoff point for determining an outlier
#' @export
#' @importFrom stats IQR quantile
upper_outlier <- function (x, extreme = 2) {
  lowerq <- quantile(x)[2]
  upperq <- quantile(x)[4]
  iqr <- IQR(x)
  x > iqr * extreme + upperq
}
