#' Plot specimens in tangent space
#'
#' This function plots a set of Procrustes-aligned specimens in tangent space along their principal axes. This is a custom replacement for `geomorph::plotTangentSpace` tailored for the fish reliability project.
#'
#' @param coords An array (p x k x n) containing landmark coordinates for a set of aligned specimens, possibly from `geomorph::gpagen`. Array labels must be in the format `group1_group2_individual`.
#' @param axis1 The name of the first axis to plot
#' @param axis2 The name of the second axis to plot
#' @return A `ggplot` object
#' @export
#' @import dplyr
#' @import ggplot2
plot_tangent_space <- function(coords, axis1="PC1", axis2="PC2") {
  x <- geomorph::two.d.array(coords)
  pc.res <- stats::prcomp(x)
  pcdata <- pc.res$x

  res <- as.data.frame(pcdata)
  res$id <- rownames(res)
  res2 <- res %>% tidyr::separate(id, into=c("family", "role","sha"), sep="_") %>% left_join(full_family, by="family")

  # consensus in Procrustes space for downstream highlighting
  dots <- list(
    x=lazyeval::interp(~median(axis1, na.rm=T), axis1=as.name(axis1)),
    y=lazyeval::interp(~median(axis2, na.rm=T), axis2=as.name(axis2)))
  consensus <- res2 %>% group_by(family, role) %>% summarise_(.dots=dots)

  chulls <- res2 %>% group_by(family) %>% do(.[chull(.[[axis1]], .[[axis2]]), ])

  ggplot(res2, aes_string(axis1, axis2)) + geom_polygon(aes_string(axis1, axis2, fill="family", color="family"), data=chulls, alpha=0.5) + geom_point(aes(color=role)) + geom_path(aes(x, y, group=family), data=consensus) + geom_point(aes(x, y), data=consensus, size=3, color="black") + geom_point(aes(x, y, color=role), data=consensus, size=2)	+ geom_dl(aes(label=family), data=res2, list("smart.grid")) + theme_minimal() + theme(legend.position="none")
}
