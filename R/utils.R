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


#' Parallel version of BAMMtools getRateThroughTimeMatrix
#'
#' @param ... other parameters, see \link[BAMMtools]{getRateThroughTimeMatrix}
#' @param mc.cores number of cores to initialize a cluster with
#' @export
par_get_rtt_matrix <- function(ephy, start.time=NULL, end.time=NULL, nslices=100, node=NULL, nodetype = 'include', mc.cores = getOption("mc.cores", 2L)) {

  if (!'bammdata' %in% class(ephy)) {
    stop("Object ephy must be of class 'bammdata'\n");
  }
  if (ephy$type == 'diversification') {
    stop("'diversification' not supported")
  }


  if (is.null(node)) {
    nodeset <- ephy$edge[,2];
  } else if (!is.null(node) & nodetype == 'include') {
    #nodeset <- getDesc(ephy, node)$desc_set;
    nodeset <- unlist(sapply(node, function(x) getDesc(ephy, x)$desc_set))
  } else if (!is.null(node) & nodetype == 'exclude') {
    nodeset <- setdiff( ephy$edge[,2], unlist(sapply(node, function(x) getDesc(ephy, x)$desc_set)));
  } else {
    stop('error in getRateThroughTimeMatrix\n');
  }

  bt <- branching.times(as.phylo.bammdata(ephy));
  maxpossible <- max(bt[as.character(intersect(nodeset, ephy$edge[,1]))]);

  if (is.null(start.time)) {
    start.time <- max(bt) - maxpossible;
  }
  if (is.null(end.time)) {
    end.time <- max(bt);
  }

  decimals = function(x){
    if(x%%1 != 0)	{
      return(nchar(strsplit(as.character(x),".",fixed=TRUE)[[1]][[2]]));
    }		else	{
      return(10);
    }
  }

  tvec <- seq(start.time, end.time, length.out= nslices);
  tol = 1*10^-decimals(ephy$eventBranchSegs[[1]][1,2]);
  mm <- matrix(NA, nrow=length(ephy$eventBranchSegs), ncol=length(tvec));

  # this is much faster than 'safeCompare'.
  goodTime <- function (vec, val, tol) {
    (vec[,2] - val <= tol) & (val - vec[,3] <= tol);
  }

  res <- mclapply(mc.cores = mc.cores, 1:length(ephy$eventBranchSegs), function(i) {
    es <- ephy$eventBranchSegs[[i]];
    events <- ephy$eventData[[i]];

    isGoodNode <- rep(TRUE, nrow(es));

    sapply(1:length(tvec), function(k) {
      isGoodTime <- goodTime(es, tvec[k], tol=tol);
      if (!(is.null(node))) { # only enter this if not the root. otherwise, only have to set once per i.
        isGoodNode <- es[,1] %in% nodeset;
      }

      estemp <- es[isGoodTime & isGoodNode, ];
      tvv <- tvec[k] - events$time[estemp[,4]];
      rates <- exponentialRate(tvv, events$lam1[estemp[,4]], events$lam2[estemp[,4]]);
      mean(rates);
    })
  })

  mm <- do.call(rbind, res)

  obj <- list();
  obj$beta <- mm;
  obj$times <- tvec;

  class(obj) <- 'bamm-ratematrix';
  obj$type = 'trait';
  return(obj);
}

