#' Do PCA on Procrustes-aligned specimens
#'
#' @param coords An array (p x k x n) containing landmark coordinates for a set of aligned specimens, possibly from \link{geomorph::gpagen}.
#' @param rename_rownames Rename and remove the rownames? If TRUE (the default), converts rownames to a column in the data frame called `id`. If a character string, instead calls that column by that string. If FALSE, preserves the rownames as-is.
#' @return an object of class `tbl_df`, `tbl`, and `data.frame`, with columns PC1 to PC...
#' @export
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
plot_tangent_space <- function(coords_df, axis1="PC1", axis2="PC2") {
  # consensus in Procrustes space for downstream highlighting
  dots <- list(
    x=lazyeval::interp(~median(axis1, na.rm=T), axis1=as.name(axis1)),
    y=lazyeval::interp(~median(axis2, na.rm=T), axis2=as.name(axis2)))
  consensus <- coords_df %>% group_by(family, role) %>% summarise_(.dots=dots)

  chulls <- coords_df %>% group_by(family) %>% do(.[chull(.[[axis1]], .[[axis2]]), ])

  ggplot(coords_df, aes_string(axis1, axis2)) + geom_polygon(aes_string(axis1, axis2, fill="family", color="family"), data=chulls, alpha=0.5) + geom_point(aes(color=role)) + geom_path(aes(x, y, group=family), data=consensus) + geom_point(aes(x, y), data=consensus, size=3, color="black") + geom_point(aes(x, y, color=role), data=consensus, size=2)	+ theme_minimal() + guides(color=FALSE)
}

#' Ungrouped plot tangent space
#'
#' XXX: Make plot tangent space take a grouped_df
#' @export
plot_tangent_space2 <- function(coords_df, axis1="PC1", axis2="PC2", groupcol="family", labelmethod=list("smart.grid")) {
  chulls <- coords_df %>% group_by_(groupcol) %>% do(.[chull(.[[axis1]], .[[axis2]]), ])

  ggplot(coords_df, aes_string(axis1, axis2, fill=groupcol, colour=groupcol)) + geom_polygon(data=chulls, alpha=0.5) + geom_point() + theme_minimal() + guides(color = FALSE)
}

#' Black theme for ggplot
#' @export
theme_black <- function(base_size = 12, base_family = "Helvetica") {
  theme(
    line =               element_line(colour = "black", size = 0.5, linetype = 1,
                                      lineend = "butt"),
    rect =               element_rect(fill = "white", colour = "black", size = 0.5, linetype = 1),
    text =               element_text(family = base_family, face = "plain",
                                      colour = "black", size = base_size,
                                      hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 0.9),
    axis.text =          element_text(size = rel(0.8), colour = "white"),
    strip.text =         element_text(size = rel(0.8), colour = "white"),

    axis.line =          element_blank(),
    axis.text.x =        element_text(vjust = 1),
    axis.text.y =        element_text(hjust = 1),
    axis.ticks =         element_line(colour = "white", size = 0.2),
    axis.title =         element_text(colour = "white"),
    axis.title.x =       element_text(vjust = 1),
    axis.title.y =       element_text(angle = 90),
    axis.ticks.length =  unit(0.3, "lines"),
    axis.ticks.margin =  unit(0.5, "lines"),

    legend.background =  element_rect(colour = NA),
    legend.margin =      unit(0.2, "cm"),
    legend.key =         element_rect(fill = "black", colour = "white"),
    legend.key.size =    unit(1.2, "lines"),
    legend.key.height =  NULL,
    legend.key.width =   NULL,
    legend.text =        element_text(size = rel(0.8), colour = "white"),
    legend.text.align =  NULL,
    legend.title =       element_text(size = rel(0.8), face = "bold", hjust = 0, colour = "white"),
    legend.title.align = NULL,
    legend.position =    "right",
    legend.direction =   "vertical",
    legend.justification = "center",
    legend.box =         NULL,

    panel.background =   element_rect(fill = "black", colour = NA),
    panel.border =       element_rect(fill = NA, colour = "white"),
    panel.grid.major =   element_line(colour = "grey20", size = 0.2),
    panel.grid.minor =   element_line(colour = "grey5", size = 0.5),
    panel.margin =       unit(0.25, "lines"),

    strip.background =   element_rect(fill = "grey30", colour = "grey10"),
    strip.text.x =       element_text(),
    strip.text.y =       element_text(angle = -90),

    plot.background =    element_rect(colour = "black", fill = "black"),
    plot.title =         element_text(size = rel(1.2)),
    plot.margin =        unit(c(1, 1, 0.5, 0.5), "lines"),

    complete = TRUE
  )
}

#' Plot BAMM color ramps
#'
#' Plots a color ramp for a BAMM phylorate plot, using the BAMM rates KDE.
#'
#' @param q BAMM plot data, from \link{BAMMtools::plot.bammdata}
#' @export
#' @author Pascal Title
plotColorRamp <- function(q) {
  y <- q$colordens[,2]
  good <- y > 0
  y <- log(y[good])
  y <- y + abs(min(y))
  x <- q$colordens[,1][good]
  cols <- q$colordens[,3][good]
  plot.new()
  plot.window(xlim = c(min(0,min(x)), max(x)), ylim = c(0, max(y)))
  segments(x, y, x, 0, lend = 2, col = cols, lwd=3)
  axis(1, signif(seq(min(0,min(x)), max(x), length.out = 5), 2), xaxs = "i", cex.axis = 1, tcl = NA, mgp = c(0, 0.25, 0))
  axis(2, round(seq(0, max(y), length.out = 3), 0), las = 1, yaxs = "i", cex.axis = 1, tcl = NA, mgp = c(0, 0.25, 0))
}


#' Plot BAMM rate through time
#' @param rtt Rate through tiem matrix
#' @param shifts list of other rate through time matrices corresponding to certain nodes
#' @param xlim vector of length 2, x limits of plot
#' @param ylim vector of length 2, y limits of plot
#' @param labels whether to plot axis labels
#' @export
plot_rtt <- function(rtt, shifts = list(), xlim = NULL, ylim = NULL, labels = TRUE) {
  vartype <- if (rtt$type == "trait") "beta" else "lambda"
  if (length(shifts) > 0) {
    allshifts <- c(list(rtt), shifts)
  } else {
    allshifts <- list(rtt)
  }
  allshifts <- lapply(allshifts, function (x) {
    x$avg <- apply(x[[vartype]], 2, median)
    x$times <- rev(x$times)
    x
  })
  xlim <- if (is.null(xlim)) c(0, max(sapply(allshifts, `[[`, "times"))) else xlim
  ylim <- if (is.null(ylim)) range(sapply(allshifts, `[[`, "avg")) else ylim
  plot.new()
  plot.window(rev(xlim), ylim, log = "y")
  axis(at = axTicks(1), side = 1, labels = axTicks(1))
  axis(at = c(1e-4, 1e-3, 1e-2, 1e-1), side = 2, las = 2, labels = labels)
  lines(allshifts[[1]]$times, allshifts[[1]]$avg)
  if (length(shifts) > 0) {
    lapply(2:length(allshifts), function (x) {
      lines(allshifts[[x]]$times - min(allshifts[[x]]$times), allshifts[[x]]$avg, col="red", lwd=1.5)
    })
  }
  return(invisible(NULL))
}

#' Plot horizontal or vertical text
#' @param str string to plot
#' @export
htext <- function(str) {
  plot.new()
  text(0.5, 0.5, str, cex=1)
}

#' @export
#' @rdname htext
vtext <- function(str) {
  plot.new()
  text(0.5, 0.5, str, cex=1, srt=90)
}

