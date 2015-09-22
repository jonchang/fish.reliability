#' Coerce to a data frame
#'
#' @param x A rectangular object, with or without rownames
#' @param ... Additional columns to include
#' @return A data frame. Row names, if present in `x`, are included as a `.rownames` column
#' @export
as_df <- function(x, ...) {
  if (!is.null(rownames(x))) {
    x <- cbind(.rownames=rownames(x), x)
  }
  x <- cbind(x, ...)
  as.data.frame(x, stringsAsFactors=F)
}

#' Prediction object for `ipred`
#' @export
ip.lda <- function(object, newdata) {
  predict(object, newdata=newdata)$class
}

#' Formats a data frame for \code{\link[geomorph]{gpagen}}
#'
#' @param df a data frame. must include columns "x" and "y"
#' @param ... columns to include as labels
#' @param .dots used to work around nonstandard evaluation
#' @return A 3D matrix
#' @export
format_for_gpagen <- function(df, ...) {
  format_for_gpagen_(df, .dots = lazyeval::lazy_dots(...))
}


#' @export
#' @rdname format_for_gpagen
format_for_gpagen_ <- function(df, ..., .dots) {
  dots <- lazyeval::all_dots(.dots, ..., all_named = TRUE)
  splut <- tidyr::gather(df, variable, val, x:y) %>% tidyr::unite_("id", names(dots))
  xtabs(val ~ mark + variable + id, data=splut)
}



#' Gets the GPA and also vectors of labels
#' @export
gpa_with_vectors <- function(df, ...) {
  fmted <- format_for_gpagen_(df, .dots = lazyeval::lazy_dots(...))
  gpaed <- gpagen(fmted, ShowPlot = FALSE)

  ald <- str_split(dimnames(fmted)$id, fixed("_")) %>%
    lapply(function(x) data.frame(t(x))) %>% bind_rows %>%
    set_colnames(lazyeval::lazy_dots(...) %>% lapply(`[[`, "expr") %>% as.character)
  list(coords = gpaed$coords, vectors = ald)
}


#' Gets replicate error
#' @export
get_replicate_error <- function (df) {
  fmted <- format_for_gpagen(df, family, role, sha1mac, sequence)
  gpaed <- gpagen(fmted, ShowPlot = FALSE)

  ald <- str_split(dimnames(fmted)$id, fixed("_")) %>%
    lapply(function(x) data.frame(t(x))) %>% bind_rows %>%
    set_colnames(c("family", "role", "individual", "sequence")) %>%
    group_by(family, individual) %>% mutate(sequence = row_number(individual))

  fam <- ald$family
  role <- ald$role
  individual <- ald$individual
  sequence <- ald$sequence

  full_model <- procD.lm(gpaed$coords ~ fam * individual * ald$sequence, iter = 1)

  ((full_model$MS[1] - full_model$MS[7])/5) / (full_model$MS[7] + ((full_model$MS[1] - full_model$MS[7])/5))
}

