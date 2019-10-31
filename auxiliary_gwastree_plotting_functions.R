#!/usr/bin/R
# -------------------------------------------
# Auxiliary functions for plotting gwas tree:
# -------------------------------------------
library(circlize)
library(dendextend)
library(ComplexHeatmap)
library(dplyr)
library(viridis)

# ------------------------
# Tree plotting functions:
# ------------------------
# Reduced circleplot FN:
plot.labels.track = function(dend){
    NLAB = length(labels(dend))
    circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
                     circos.text(1:NLAB-0.5, rep(0, NLAB), labels(dend), col = labels_colors(dend), cex=.25,
                                 facing = "clockwise", niceFacing=TRUE, adj = c(0, 0.5)) }, bg.border = NA, track.height = 0.1)
}


plot.metadata.track = function(colmat, track.height=0.075, reduced.metadata=FALSE){
    if (reduced.metadata) {
        colmat = as.matrix(colmat[,1])
    }
    ncat = ncol(colmat)
    track.height = track.height / 5 * ncat
    circos.track(ylim = c(0, ncat), bg.border = NA, panel.fun = function(x, y) {
                     m = t(colmat[,ncat:1])
                     nr = nrow(m)
                     nc = ncol(m)
                     for(i in 1:nr) {
                         circos.rect(1:nc - 1, rep(nr - i, nc), 
                                     1:nc, rep(nr - i + 1, nc), 
                                     border = NA, col = m[i, ])
                     } }, track.height=track.height) 
}



plot.circos.tree = function(dend, ntdf=NULL, track.height=0.5, curved.lines=FALSE, boxdf=NULL, boxlabdf=NULL, boxhoriz=FALSE) {
    max_height = max(attr(dend, "height"))
    circos.track(ylim = c(0, max_height), panel.fun = function(x, y) {
                     # circos.dendrogram(dend, max_height = max_height)
                     circos.dendrogram.internal(dend, max_height = max_height, curved.lines=curved.lines)
                     if (!is.null(ntdf)){
                         # Separate nodes that are on the left vs. right by adj.
                         if (is.null(ntdf$cex)){ ntdf$cex = .35 }
                         ntdf$face = 0
                         ntdf$face[ntdf$X1 <= 209] = 1
                         ntdf$face[ntdf$X1 >= 626] = 1
                         ntdf$top = sapply(ntdf$X1, function(x){
                                                 abs(min(x  %% 833, (x + 418) %% 833) - 209)/209})
                         ntdf$top = ntdf$top * (2 * (ntdf$X1 < 418) - 1)
                         # diff to max height:
                         ntdf$maxy = max_height - ntdf$X2 - .05 # tolerance
                         ntdf$maxy = sapply(ntdf$maxy, function(x){min(0.5, x)})
                         ntdf$ry = runif(nrow(ntdf), min=min(.05, ntdf$maxy), 
                                         max=ntdf$maxy)
                         lndf = ntdf[ntdf$face == 0,]
                         rndf = ntdf[ntdf$face == 1,]
                         if (nrow(rndf) > 0){
                             circos.text(rndf$X1 + ux(5 * rndf$top, "mm"), (max_height - rndf$X2) - rndf$ry, rndf$symbol, 
                                         col=rndf$color, cex=rndf$cex, facing = "downward", niceFacing=TRUE, adj = c(1, 1)) 
                             # Draw line between node and text:
                             circos.segments(x0=rndf$X1, y0=max_height - rndf$X2, x1=rndf$X1 + ux(5 * rndf$top, "mm"), 
                                             y1=max_height - rndf$X2 - rndf$ry + 0.01, straight=TRUE, col=rndf$color)
                         }
                         if (nrow(lndf) > 0){
                             circos.text(lndf$X1 + ux(-5 * lndf$top, "mm"), (max_height - lndf$X2) - lndf$ry, lndf$symbol, 
                                         col=lndf$color, cex=lndf$cex, facing = "downward", niceFacing=TRUE, adj = c(0, 1)) 
                             # Draw line between node and text:
                             circos.segments(x0=lndf$X1, y0=max_height - lndf$X2, x1=lndf$X1 + ux(-5 * lndf$top, "mm"), 
                                             y1=max_height - lndf$X2 - lndf$ry + 0.01, straight=TRUE, col=lndf$color)
                         }
                     }
                     if (!is.null(boxdf)){ add.circos.boxes(boxdf, boxlabdf, boxhoriz=boxhoriz) }
                     }, track.height=track.height, bg.border = NA)
}



# Auxiliary for circosplot:
plot.rna.avail = function(rna.avail, scale=1){
    circos.track(ylim = c(0, 1), bg.border = NA, panel.fun = function(x, y) {
                     for(i in which(rna.avail == 1)) {
                         circos.rect(i - 1, 0, i, 1, 
                                     col='grey25', border=NA)
                     } }, track.height=0.0075 * scale) 
}


# Auxiliary for circosplot:
plot.hit.track = function(hit.track, scale=1){
    hpal = rev(magma(100))[10:100]
    top.nhit = max(max(hit.track$nhit), 1)
    hbins <- cut(1:top.nhit, seq(0, top.nhit, length.out=length(hpal)), include.lowest=T) 
    hcol = hpal[hbins]
    circos.track(ylim = c(0, 1), bg.border = NA, panel.fun = function(x, y) {
                     for(j in 1:nrow(hit.track)) {
                         i = hit.track$i[j]
                         nhit = hit.track$nhit[j]
                         circos.rect(i - 1, 0, i, 1, 
                                     col=hcol[nhit], border=NA)
                     } }, track.height=0.005 * scale) 
}


plot.leafpval.track = function(leaf.pval, cutp=CUTP, scale=1){
    ind = which(leaf.pval > 0)
    if (length(ind) > 0){
        lp = leaf.pval[ind]
        palette = col4
        # Color for strength of p-value:
        bins <- cut(lp, seq(0, cutp, length.out=length(palette)), include.lowest=T) 
        lcol = palette[bins]
        circos.track(ylim = c(0, 1), bg.border = NA, panel.fun = function(x, y) {
                         for(j in 1:length(ind)) {
                             i = ind[j]
                             circos.rect(i - 1, 0, i, 1, 
                                         col=lcol[j], border=NA)
                         } }, track.height=0.005 * scale) 
    }
}

plot.fractional.tracks = function(udf, fractional=TRUE, fulltrack=TRUE, scale=1){
    # Add barplots:
    nr = nrow(udf)
    udf$col = as.character(udf$col)
    if (fulltrack){
        topy = max(udf$total)
        circos.track(ylim = c(0, topy), bg.border = NA, panel.fun = function(x, y) {
                         for(i in 1:nr) {
                             circos.rect(rep(i - 1, 2),c(0, udf$uq[i]),
                                         rep(i, 2),c(udf$uq[i], udf$total[i]),
                                         col=c(udf$col[i], 'grey85'), border=NA)
                         } }, track.height=0.06 * scale) 
    }
    if (fractional){
        circos.track(ylim = c(0, 1), bg.border = NA, panel.fun = function(x, y) {
                         for(i in 1:nr) {
                             circos.rect(rep(i - 1, 2),c(0, udf$frac[i]),
                                         rep(i, 2),c(udf$frac[i], 1),
                                         col=c(udf$col[i], 'grey85'), border=NA)
                         } }, track.height=0.06 * scale) 
    } 
}


add.circos.boxes = function(boxdf, boxlabdf=NULL, boxhoriz=FALSE){
    circos.rect(boxdf$X1, boxdf$Y1,
                boxdf$X2, boxdf$Y2 - max_height * 0.005,
                lwd=boxdf$lwd, lty='dashed', col=NA,border=boxdf$col)
    if(!is.null(boxlabdf)){
        if (boxhoriz) { 
            facing='downward' 
            fix.adj=TRUE
            adj = 0.925
        } else {
            facing='reverse.clockwise'
            adj = 0.975
            fix.adj=FALSE
        }
        circos.segments(boxlabdf$X1, adj * max_height - boxlabdf$X2,
                        boxlabdf$X1, max_height - boxlabdf$X2)
        par(xpd=NA)
        if (fix.adj){
            boxlabdf$hadj = (cos((NLAB - boxlabdf$X1)/NLAB * 2* pi) + 1 ) / 2
            boxlabdf$vadj = (sin((NLAB - boxlabdf$X1)/NLAB * 2* pi) + 1 ) / 2
            for (i in 1:nrow(boxlabdf)) {
                circos.text(boxlabdf$X1[i], adj * max_height - boxlabdf$X2[i], as.character(boxlabdf$node)[i], 
                            col=boxdf$col[i], cex=.75, facing=facing, niceFacing=TRUE, adj=c(boxlabdf$hadj[i], boxlabdf$vadj[i])) 
            }
        } else {
            circos.text(boxlabdf$X1, adj * max_height - boxlabdf$X2, as.character(boxlabdf$node), 
                        col=boxdf$col, cex=.75, facing=facing, niceFacing=TRUE, adj=c(0, 0.5)) 
        }
        par(xpd=FALSE)
    } 
}

# Full plotting function:
circleplot = function(dend, lab, with.tree=TRUE, nodedf=NULL, udf=NULL,
                      # For plotting metadata:
                      add.metadata=TRUE, hit.track=NULL, scale=1,
                      reduced.metadata=FALSE, curved.lines=FALSE,
                      leaf.pval=NULL, rna.avail=NULL, plot.labels=TRUE,
                      # For plotting tracks around:
                      fractional=TRUE, fulltrack=TRUE, cutp=CUTP, 
                      boxdf=NULL, boxlabdf=NULL, boxhoriz=FALSE){
    # Set up dendrogram:
    NLAB = length(labels(dend))
    mll = meta.image(metamat[lab,5:1], colvals=colvals, cex=.6, horiz=TRUE, return.mat=TRUE)
    nummat = mll[[1]]
    colsmeta = mll[[2]]
    colmat = matrix(colsmeta[nummat], nrow=nrow(nummat))
    labels_dend <- labels(dend)
    if (as.logical(anyDuplicated(labels_dend))) {
        labels(dend) <- paste0(seq_along(labels_dend), "_", labels_dend)
        labels_dend <- labels(dend) }
    # Initialize figure:
    circos.par(cell.padding = c(0, 0, 0, 0), track.margin=rep(0.001,2))
    circos.initialize(factors = "single", xlim = c(0, NLAB)) 
    # Plot tracks:
    if (plot.labels){ plot.labels.track(dend) }
    if (!is.null(hit.track)){ plot.hit.track(hit.track, scale=scale) } 
    if (!is.null(leaf.pval)){ plot.leafpval.track(leaf.pval, cutp, scale=scale) } 
    if (!is.null(udf)){ plot.fractional.tracks(udf, fractional=fractional, fulltrack=fulltrack, scale=1) }
    if (!is.null(rna.avail)){ plot.rna.avail(rna.avail, scale=scale) } 
    if (add.metadata){ plot.metadata.track(colmat, track.height=0.075 * scale, reduced.metadata=reduced.metadata) }
    if (with.tree){ plot.circos.tree(dend, ntdf=nodedf, track.height=0.5, curved.lines=curved.lines, boxdf=boxdf, boxlabdf=boxlabdf, boxhoriz=boxhoriz) }
    circos.clear()
    circos.clear()
}




# Same as above, but split:
circleplot_split = function(dend, lab, NCLUST, with.tree=TRUE, only.diff=FALSE, nodedf=NULL, curved.lines=FALSE){
    # Set up dendrogram:
    NLAB = length(labels(dend))
    mll = meta.image(metamat[lab,5:1], colvals=colvals, cex=.6, horiz=TRUE, return.mat=TRUE)
    nummat = mll[[1]]
    colsmeta = mll[[2]]
    colmat = matrix(colsmeta[nummat], nrow=nrow(nummat))
    labels_dend <- labels(dend)
    if (as.logical(anyDuplicated(labels_dend))) {
        labels(dend) <- paste0(seq_along(labels_dend), "_", labels_dend)
        labels_dend <- labels(dend) }
    # Plot figure:
    # Make split dendrogram:
    colpair = colorRampPalette(brewer.pal(n=12,name="Paired"))(NCLUST)
    dend <- color_branches(dend, k=NCLUST, col=colpair)
    dend = set(dend, "labels_cex", .18)
    labels_dend <- labels(dend)
    if (as.logical(anyDuplicated(labels_dend))) {
        labels(dend) <- paste0(seq_along(labels_dend), "_", labels_dend)
        labels_dend <- labels(dend)
    }
    # Factors/splits:
    kcol = unlist(get_leaves_attr(dend, 'edgePar'))
    letext = c(letters, paste0('a', letters))
    factors = letext[as.numeric(factor(kcol, levels = unique(kcol)))]
    dend_list = get_subdendrograms(dend, k=NCLUST)
    # Correct the dendrogram order:
    dnames = sapply(dend_list, function(x){as.numeric(strsplit(labels(x)[1],'_')[[1]][[1]]) })
    names(dend_list) = letext[1:NCLUST][order(order(dnames))]
    # ----------------
    # Plot dendrogram:
    circos.par(cell.padding = c(0, 0, 0, 0), track.margin=rep(0.001,2))
    circos.initialize(factors, xlim = cbind(c(0, 0), table(factors)))
    # Text labels:
    circos.track(ylim = c(0, 1), 
                 panel.fun = function(x, y) {
                     sector.index = CELL_META$sector.index
                     dend = dend_list[[sector.index]]
                     txt = labels(dend)
                     nc = length(txt)
                     circos.text(1:nc-0.5, rep(0, nc), txt, 
                                 col = labels_colors(dend), cex=.25,
                                 facing = "clockwise", niceFacing=TRUE, adj = c(0, 0.5)) 
                 }, bg.border = NA, track.height = 0.1)
    # Metadata:
    circos.track(ylim = c(0, 5), bg.border = NA, panel.fun = function(x, y) {
                     sector.index = CELL_META$sector.index
                     m = t(colmat[factors == sector.index,5:1])
                     nr = nrow(m)
                     nc = ncol(m)
                     for(i in 1:nr) {
                         circos.rect(1:nc - 1, rep(nr - i, nc), 
                                     1:nc, rep(nr - i + 1, nc), 
                                     border = NA, col = m[i, ])
                 } }, track.height=0.125)
    # Dendrogram
    max_height = max(sapply(dend_list, function(x) attr(x, "height")))
    circos.track(ylim = c(0, max_height), bg.border = NA, track.height = 0.3, 
                 panel.fun = function(x, y) {
                     sector.index = CELL_META$sector.index
                     dend = dend_list[[sector.index]]
                     # circos.dendrogram(dend, max_height = max_height)
                     circos.dendrogram.internal(dend, max_height = max_height, curved.lines=curved.lines)
                     if (!is.null(nodedf)){
                         circos.text(nodedf$V1, max_height - nodedf$V2, as.character(nodedf$node),
                                     col = 'black', cex=.25,
                                     facing = "reverse.clockwise", niceFacing=TRUE, adj = c(0, 0.5)) 
                     }
                 }
    )
}


calculate.curve = function(x1,x2, y1,y2, curve=NULL) {
    if (is.null(curve)){
        # Curve for connecting line:
        sq = seq(-4,max(y1 - y2 * 20, 8), length.out=50)
        curve = 1/ (1 + exp(-sq))
        curve = (curve - min(curve)) / (diff(range(curve)))
    }
    x = (curve * (x1 - x2)) + x2
    y = seq(y2, y1, length.out=length(curve))
    return(list(x=x, y=y))
}

# ------------------------------------------------------------------------------
# Dendrogram plotting:
# Same as function in: https://github.com/jokergoo/circlize/blob/master/R/plot.R
# Changes: Skips plotting lty = 0 (no-line)
# ------------------------------------------------------------------------------
circos.dendrogram.internal = function(dend, facing = c("outside", "inside"), 
                                      max_height = NULL, use_x_attr = FALSE, 
                                      curved.lines=FALSE) {
	facing = match.arg(facing)[1]
    if(is.null(max_height)) { max_height = attr(dend, "height") }
    is.leaf = function(object) {
        leaf = attr(object, "leaf")
        if(is.null(leaf)) { FALSE } else { leaf }
    }

    use_x_attr = use_x_attr
    lines_par = function(col = par("col"), lty = par("lty"), lwd = par("lwd"), ...) { return(list(col = col, lty = lty, lwd = lwd)) }
    points_par = function(col = par("col"), pch = par("pch"), cex = par("cex"), ...) { return(list(col = col, pch = pch, cex = cex)) }

    draw.d = function(dend, max_height, facing = "outside", max_width = 0) {
        leaf = attr(dend, "leaf")
        d1 = dend[[1]]  # child tree 1
        d2 = dend[[2]]  # child tree 2
        height = attr(dend, "height")
        midpoint = attr(dend, "midpoint")

        if(use_x_attr) {
            x1 = attr(d1, "x")	
        } else {
            if(is.leaf(d1)) {
                x1 = x[as.character(attr(d1, "label"))]
            } else {
                x1 = attr(d1, "midpoint") + x[as.character(labels(d1))[1]]
            }
        }
        y1 = attr(d1, "height")

        if(use_x_attr) {
            x2 = attr(d2, "x")	
        } else {
            if(is.leaf(d2)) {
                x2 = x[as.character(attr(d2, "label"))]
            } else {
                x2 = attr(d2, "midpoint") + x[as.character(labels(d2))[1]]
            }
        }
        y2 = attr(d2, "height")

        # graphic parameter for current branch
        # only for lines, there are lwd, col, lty
        edge_par1 = do.call("lines_par", as.list(attr(d1, "edgePar")))  # as.list to convert NULL to list()
        edge_par2 = do.call("lines_par", as.list(attr(d2, "edgePar")))
        node_par = attr(dend, "nodePar")
        if(!is.null(node_par)) node_par = do.call("points_par", as.list(attr(dend, "nodePar")))

        if (curved.lines){
            if(facing == "outside") {
                if (edge_par1$lty != 0){
                    lxy = calculate.curve(x1=x1, x2=(x1 + x2)/2, y1=max_height - y1,
                                           y2=max_height - height, curve=NULL)
                    circos.lines(lxy$x, lxy$y, col=edge_par1$col, lty=edge_par1$lty, 
                                 lwd=edge_par1$lwd, straight=TRUE)
                                 # lwd=0.1, straight=TRUE)
                }
                if (edge_par2$lty != 0){
                    lxy = calculate.curve(x1=x2, x2=(x1 + x2)/2, y1=max_height - y2,
                                           y2=max_height - height, curve=NULL)
                    circos.lines(lxy$x, lxy$y, col=edge_par2$col, lty=edge_par2$lty,
                                 lwd=edge_par2$lwd, straight=TRUE)
                                 # lwd=0.1, straight=TRUE)
                }
                if(!is.null(node_par)) {
                    circos.points((x1+x2)/2, max_height - height, col = node_par$col, 
                                  pch = node_par$pch, cex = node_par$cex)
                }
            } else if(facing == "inside") {
                if (edge_par1$lty != 0){
                    lxy = calculate.curve(x1=x1, x2=(x1 + x2)/2, 
                                           y1=y1, y2=height, curve=NULL)
                    circos.lines(lxy$x, lxy$y, col=edge_par1$col, lty=edge_par1$lty, 
                                 lwd=edge_par1$lwd, straight=TRUE)
                }
                if (edge_par2$lty != 0){
                    lxy = calculate.curve(x1=x2, x2=(x1 + x2)/2, 
                                           y1=y2, y2=height, curve=NULL)
                    circos.lines(lxy$x, lxy$y, col=edge_par2$col, lty=edge_par2$lty,
                                 lwd=edge_par2$lwd, straight=TRUE)
                }
                if(!is.null(node_par)) {
                    circos.points((x1+x2)/2, height, col = node_par$col,
                                  pch = node_par$pch, cex = node_par$cex)
                }
            }
        } else {
            if(facing == "outside") {
                if (edge_par1$lty != 0){
                    circos.lines(c(x1, x1), max_height - c(y1, height), col=edge_par1$col, 
                                 lty=edge_par1$lty, lwd=edge_par1$lwd, straight=TRUE)
                    circos.lines(c(x1, (x1+x2)/2), max_height - c(height, height), col=edge_par1$col, 
                                 lty=edge_par1$lty, lwd=edge_par1$lwd)
                }
                if (edge_par2$lty != 0){
                    circos.lines(c(x2, x2), max_height - c(y2, height), col=edge_par2$col, 
                                 lty=edge_par2$lty, lwd=edge_par2$lwd, straight=TRUE)
                    circos.lines(c(x2, (x1+x2)/2), max_height - c(height, height), col=edge_par2$col,
                                 lty=edge_par2$lty, lwd=edge_par2$lwd)
                }
                if(!is.null(node_par)) {
                    circos.points((x1+x2)/2, max_height - height, col = node_par$col, 
                                  pch = node_par$pch, cex = node_par$cex)
                }
            } else if(facing == "inside") {
                if (edge_par1$lty != 0){
                    circos.lines(c(x1, x1), c(y1, height), col = edge_par1$col,
                                 lty = edge_par1$lty, lwd = edge_par1$lwd, straight = TRUE)
                    circos.lines(c(x1, (x1+x2)/2), c(height, height), col = edge_par1$col,
                                 lty = edge_par1$lty, lwd = edge_par1$lwd)
                }
                if (edge_par2$lty != 0){
                    circos.lines(c(x2, x2), c(y2, height), col = edge_par2$col,
                                 lty = edge_par2$lty, lwd = edge_par2$lwd, straight = TRUE)
                    circos.lines(c(x2, (x1+x2)/2), c(height, height), col = edge_par2$col,
                                 lty = edge_par2$lty, lwd = edge_par2$lwd)
                }
                if(!is.null(node_par)) {
                    circos.points((x1+x2)/2, height, col = node_par$col,
                                  pch = node_par$pch, cex = node_par$cex)
                }
            }
        }
        # do it recursively
        if(is.leaf(d1)) {
            node_par = attr(d1, "nodePar")
            if(!is.null(node_par)) node_par = do.call("points_par", as.list(attr(d1, "nodePar")))
            if(facing == "outside") {
                if(!is.null(node_par)) {
                    circos.points(x1, max_height, col = node_par$col, pch = node_par$pch, cex = node_par$cex)
                }
            } else if(facing == "inside") {
                if(!is.null(node_par)) {
                    circos.points(x1, 0, col = node_par$col, pch = node_par$pch, cex = node_par$cex)
                }
            }
        } else {
            draw.d(d1, max_height, facing, max_width)
        }

        if(is.leaf(d2)) {
            node_par = attr(d2, "nodePar")
            if(!is.null(node_par)) node_par = do.call("points_par", as.list(attr(d2, "nodePar")))
            if(facing == "outside") {
                if(!is.null(node_par)) {
                    circos.points(x2, max_height, col = node_par$col, pch = node_par$pch, cex = node_par$cex)
                }
            } else if(facing == "inside") {
                if(!is.null(node_par)) {
                    circos.points(x2, 0, col = node_par$col, pch = node_par$pch, cex = node_par$cex)
                }
            }
        } else {
            draw.d(d2, max_height, facing, max_width)
        }
    }
    labels = as.character(labels(dend))
    x = seq_along(labels) - 0.5
    names(x) = labels
    n = length(labels)
    draw.d(dend, max_height, facing, max_width = n)
}
