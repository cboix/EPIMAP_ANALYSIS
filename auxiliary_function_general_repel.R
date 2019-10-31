#!/usr/bin/R
# -----------------------------------
# Auxiliary function for label repel:
# -----------------------------------
# library(ggrepel)
# library(rlang)
# importFrom(grid,convertHeight)
# importFrom(grid,convertWidth)

#' Repel labels from points and each other
#'
#' @param x
#' @param y
#' @param labels
#' @return Return new locations and facings for each label
#' @noRd
general_repel_text = function(x, y, labels, cex=1, pch=19, pt.cex=1, 
                              seed=NULL, hjust=0.5, vjust=0.5,
                              force = 1, force_pull = 1, max.iter = 2000,
                              xlim = c(NA, NA), ylim = c(NA, NA), direction = "both"){
    require(ggrepel)
    require(rlang)
    if (length(hjust) != length(x)){ hjust = rep_len(hjust, length.out=length(x)) }
    if (length(vjust) != length(x)){ vjust = rep_len(vjust, length.out=length(x)) }
    # --------------
    # Label padding:
    # --------------
    # Get strheight/strwidth
    box_widths = strwidth(labels, cex=cex)
    box_heights = strheight(labels, cex=cex)
    # The padding around each bounding box (for now, based on the height)
    box_padding_x <- .1 * mean(box_heights)
    box_padding_y <- .1 * mean(box_heights)
    # box_padding_x <- convertWidth(x$box.padding, "native", valueOnly = TRUE)
    # box_padding_y <- convertHeight(x$box.padding, "native", valueOnly = TRUE)

    # --------------
    # Point padding:
    # --------------
    # if (is.na(x$point.padding)) { x$point.padding = unit(0, "lines") }
    # point_padding_x <- convertWidth(x$point.padding, "native", valueOnly = TRUE)
    # point_padding_y <- convertHeight(x$point.padding, "native", valueOnly = TRUE)
    point_padding_x <- strheight('a', cex=pt.cex)
    point_padding_y <- strheight('a', cex=pt.cex)
    # The padding around each point.
    # if (length(x$point.size) == 1 && is.na(x$point.size)) { x$point.size = unit(0, "lines") }
    # point_size <- convertWidth(x$point.size, "native", valueOnly = TRUE)
    # point_size <- convertWidth(x$point.size, "native", valueOnly = TRUE)
    point_size = strheight('a', cex=pt.cex)

    # Do not create text labels for empty strings.
    valid_strings <- which(ggrepel:::not_empty(labels))
    invalid_strings <- which(!ggrepel:::not_empty(labels))

    # Create a dataframe with x1 y1 x2 y2
    boxes <- lapply(valid_strings, function(i) {
                        hj <- hjust[i]
                        vj <- vjust[i]
                        gw = box_widths[i]
                        gh = box_heights[i]
                        c("x1" = x[i] - gw * hj - box_padding_x,
                          "y1" = y[i] - gh * vj - box_padding_y,
                          "x2" = x[i] + gw * (1 - hj) + box_padding_x,
                          "y2" = y[i] + gh * (1 - vj) + box_padding_y) })

    # Make the repulsion reproducible if desired.
    if (is.null(seed) || !is.na(seed)) { set.seed(seed) }
    points_valid_first <- cbind(c(x[valid_strings], x[invalid_strings]),
                                c(y[valid_strings], y[invalid_strings]))
    if (length(point_size) != length(x)) {
        point_size <- rep_len(point_size, length.out=length(x)) 
    }
    point_size <- c(point_size[valid_strings], point_size[invalid_strings])

    # Repel overlapping bounding boxes away from each other.
    # repel <- ggrepel:::repel_boxes2(
    # .Call("_ggrepel_repel_boxes", PACKAGE = "ggrepel", data_point
    #       point_padding_x, point_padding_y, boxes, xlim, ylim,
    #       hjust, vjust, force, maxiter, direction)
    repel <- ggrepel:::repel_boxes(data_points = points_valid_first,
                                   # point_size = point_size,
                                   point_padding_x = point_padding_x,
                                   point_padding_y = point_padding_y,
                                   boxes = do.call(rbind, boxes),
                                   # Params:
                                   xlim = range(xlim),
                                   ylim = range(ylim),
                                   hjust = hjust %||% 0.5,
                                   vjust = vjust %||% 0.5,
                                   # NOTE: Defaults from ggrepel
                                   force = force * 1e-6,
                                   # force_pull = force_pull * 1e-2,
                                   maxiter = max.iter,
                                   direction = direction)
    return(data.frame(x=repel$x, y=repel$y,
                      # Position of original data points.
                      x.orig=x[valid_strings], y.orig=y[valid_strings], 
                      lab=labels[valid_strings]))
}
