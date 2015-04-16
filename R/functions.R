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

#' Formats a data frame for `geomorph::gpagen`
#'
#' @param df a data frame, with certain columns
#' @return A 3D matrix
#' @importFrom tidyr gather unite
#' @importFrom magrittr %>%
#' @export
format_for_gpagen <- function(df) {
  splut <- gather(df, variable, val, x:y) %>% unite(col=id, Family, role, sha1mac)
  xtabs(val ~ mark + variable + id, data=splut)
}


#' Hack the shit out of geomorph's slow procD.lm function
#' @export
my.procD.lm <- function (f1, iter = 999, RRPP = FALSE, int.first = FALSE, verbose = FALSE)
{
  form.in <- formula(f1)
  Y <- eval(form.in[[2]], parent.frame())
  if (length(dim(Y)) == 3)
    Y <- two.d.array(Y)
  else Y <- as.matrix(Y)
  if (nrow(Y) != nrow(model.frame(form.in[-2])))
    stop("Different numbers of specimens in dependent and independent variables")
  form.new <- as.formula(paste(c("Y", form.in[[3]]), collapse = "~"))
  if (int.first == TRUE)
    ko = TRUE
  else ko = FALSE
  Terms <- terms(form.new, keep.order = ko)
  mod.mf <- model.frame(Terms)
  if (any(is.na(Y)) == T)
    stop("Response data matrix (shape) contains missing values. Estimate these first (see 'estimate.missing').")
  anova.parts.obs <- geomorph:::anova.parts(f1 = form.new, Yalt = "observed",
                                 keep.order = ko)
  anova.tab <- anova.parts.obs$table
  Xs <- geomorph:::mod.mats(form.new, mod.mf, keep.order = ko)
  k <- length(Xs$Xs) - 1
  P <- array(0, c(k, 1, iter + 1))
  SS.obs <- anova.parts.obs$SS[1:k]
  P[, , 1] <- SS.obs
  method <- ifelse(RRPP, "RRPP", "resample")
  result <- simplify2array(parallel::mclapply(1:iter, function(x) geomorph:::SS.random(Y, Xs, SS.obs, Yalt=method)$SS, mc.cores=parallel::detectCores()))
  P[, , 1:iter+1] <- result
  P.val <- geomorph:::Pval.matrix(P)
  Z <- geomorph:::Effect.size.matrix(P)
  anova.tab <- data.frame(anova.tab, Z = c(Z, NA, NA), P.value = c(P.val,
                                                                   NA, NA))
  if (RRPP == TRUE) {
    anova.title = "\nRandomized Residual Permutation Procedure used\n"
  }
  else anova.title = "\nRandomization of Raw Values used\n"
  attr(anova.tab, "heading") <- paste("\nType I (Sequential) Sums of Squares and Cross-products\n",
                                      anova.title)
  class(anova.tab) <- c("anova", class(anova.tab))
  if (verbose == TRUE) {
    list(anova.table = anova.tab, call = match.call(), SS.rand = P)
  }
  else anova.tab
}
