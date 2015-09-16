#' Do template substitution on a file
#'
#' Given a file with strings that look like \code{\{templ\}}, replace it
#'
#' @param template string with templates that look like \code{\{template\}}
#' @param ... substitutions to perform
#' @param .dots Used to work around non-standard evaluation. See
#' \code{vignette("nse", "dplyr")} for details
#' @return a character vector
#' @export
#' @examples
#' substitute_template("{animal} says {noise}", animal="dog", noise="bark")
substitute_template <- function(template, ...) {
  substitute_template_(template, .dots = lazyeval::lazy_dots(...))
}

#' @export
#' @rdname substitute_template
substitute_template_ <- function(template, ..., .dots) {
  dots <- lazyeval::all_dots(.dots, ..., all_named = TRUE)
  for (name in names(dots)) {
    value <- lazyeval::lazy_eval(dots[[name]])
    template <- str_replace_all(template, fixed(paste0("{", name, "}")), value)
  }
  template
}


#' Calculates sampling fractions for BAMM
#'
#' @param phy an object of class \code{"phylo"}
#' @param grouping the grouping level at which to calculate sampling fraction
#' @param richness a data frame, with columns of groups (e.g., family, genus)
#' and a final column \code{richness} specifying how many species are in that group
#' @param tips that should be retained in the tree, but ignored for the
#' purposes of calculating taxon sampling.
#' @return a list with two components:
#' \describe{
#'   \item{global}{a numeric vector of length 1; the global sampling fraction
#'   in \code{phy}}
#'   \item{clade}{a data frame with three columns: tip labels, the grouping,
#'   and the sampling fraction}
#' }
#' @seealso \code{\link{make_sampling_file}} to write a sampling file in BAMM's format
#' @export
calculate_sampling <- function(phy, grouping, richness, ignore.tips=c()) {
  calculate_sampling_(phy, lazyeval::lazy(grouping), richness, ignore.tips=c())
}

#' @rdname calculate_sampling
#' @export
calculate_sampling_ <- function(phy, grouping, richness, ignore.tips=c()) {
  quot <- as.character(grouping$expr)
  tips <- data_frame(tip = phy$tip.label) %>%
    separate(tip, c("genus", "epithet"), sep="_", remove = FALSE, extra = "merge")

  n_tips <- nrow(tips)

  global_frac <- tips %>% distinct_(grouping) %>%
    left_join(richness, by=quot) %>%
    summarise(frac = min(n_tips / sum(richness), 1))

  clade_richness <- tips %>% left_join(richness, by=quot) %>%
    group_by_(grouping) %>% mutate(frac = min(n() / richness, 1))

  list(global = global_frac[[1]],
       clade = clade_richness %>% select_(~tip, grouping, ~frac))
}

#' Writes sampling fractions for BAMM
#'
#' @param sampling a list with two components, from \code{\link{make_sampling}}
#' @param file  A \link[base]{connection}, or a character string naming the file to print to. If "" (the default), prints to the standard output connection, the console unless redirected by \code{\link[base]{sink}}.
#' @return None (invisible \code{NULL})
#' @export
make_sampling_file <- function(sampling, file="") {
  if (is.character(file)) {
    if (file == "") file <- stdout()
  } else {
    file <- file(file, "w")
    on.exit(close(file))
  }
  cat(sampling$global, file = file, fill = TRUE)
  write.table(sampling$clade, file, append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  invisible(NULL)
}
