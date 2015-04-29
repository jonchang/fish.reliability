#' Do PCA on Procrustes-aligned specimens
#'
#' @param coords An array (p x k x n) containing landmark coordinates for a set of aligned specimens, possibly from \link{geomorph::gpagen}.
#' @param rename_rownames Rename and remove the rownames? If TRUE (the default), converts rownames to a column in the data frame called `id`. If a character string, instead calls that column by that string. If FALSE, preserves the rownames as-is.
#' @return an object of class `tbl_df`, `tbl`, and `data.frame`, with columns PC1 to PC...
#' @export
#' @import dplyr
pca_df <- function(coords, rename_rownames = TRUE) {
  x <- geomorph::two.d.array(coords)
  pc.res <- stats::prcomp(x)
  pcdata <- pc.res$x

  res <- as.data.frame(pcdata)
  if (rename_rownames) {
    if (rename_rownames == TRUE) rename_rownames <- "id"
    res[[rename_rownames]] <- rownames(res)
    rownames(res) <- NULL
  }
  tbl_df(res)
}

#' Plot specimens in tangent space
#'
#' This function plots a set of Procrustes-aligned specimens in tangent space along their principal axes. This is a custom replacement for \code{geomorph::plotTangentSpace} tailored for the fish reliability project.
#'
#' @param coords_df A PCA data frame, possibly from \link{pca_df}. Array labels must be in the format \code{group1_group2_individual}.
#' @param axis1 The name of the first axis to plot
#' @param axis2 The name of the second axis to plot
#' @return A \code{ggplot} object
#' @export
#' @import dplyr
#' @import ggplot2
plot_tangent_space <- function(coords_df, axis1="PC1", axis2="PC2") {
  res2 <- coords_df %>% tidyr::separate(id, into=c("family", "role","sha"), sep="_")

  # consensus in Procrustes space for downstream highlighting
  dots <- list(
    x=lazyeval::interp(~median(axis1, na.rm=T), axis1=as.name(axis1)),
    y=lazyeval::interp(~median(axis2, na.rm=T), axis2=as.name(axis2)))
  consensus <- res2 %>% group_by(family, role) %>% summarise_(.dots=dots)

  chulls <- res2 %>% group_by(family) %>% do(.[chull(.[[axis1]], .[[axis2]]), ])

  ggplot(res2, aes_string(axis1, axis2)) + geom_polygon(aes_string(axis1, axis2, fill="family", color="family"), data=chulls, alpha=0.5) + geom_point(aes(color=role)) + geom_path(aes(x, y, group=family), data=consensus) + geom_point(aes(x, y), data=consensus, size=3, color="black") + geom_point(aes(x, y, color=role), data=consensus, size=2)	+ geom_dl(aes(label=family), data=res2, list("smart.grid")) + theme_minimal() + theme(legend.position="none")
}

#' Ungrouped plot tangent space
#'
#' XXX: Make plot tangent space take a grouped_df
#' @export
#' @import dplyr
#' @import grDevices
#' @import directlabels
plot_tangent_space2 <- function(coords_df, axis1="PC1", axis2="PC2", groupcol="family", labelmethod=list("smart.grid")) {
  chulls <- coords_df %>% group_by_(groupcol) %>% do(.[chull(.[[axis1]], .[[axis2]]), ])

  ggplot(coords_df, aes_string(axis1, axis2, fill=groupcol, colour=groupcol)) + geom_polygon(data=chulls, alpha=0.5) + geom_point() + geom_dl(aes_string(label=groupcol), labelmethod) + theme_minimal()
}

