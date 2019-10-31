#!/usr/bin/R
# -------------------------------------
# Auxiliary functions for epimap paper:
# -------------------------------------

#' Get percentiles for each value
#'
#' @param x
#' @return percentile rank
#' @noRd
perc.rank <- function(x) trunc(rank(x))/length(x)

#' Clamp a matrix at certain quantiles
#'
#' @param mat
#' @param alpha optional cutoff.
#' @return clamped matrix
#' @noRd
clamp.mat <- function(mat, alpha=0.01){
    quant = quantile(mat, c(0.01, 1 - 0.01), na.rm=T)
    mat[mat < quant[1]] = quant[1]
    mat[mat > quant[2]] = quant[2]
    return(mat)
}

#' Clamp a matrix at certain quantiles
#'
#' @param mat
#' @param alpha optional cutoff.
#' @return clamped matrix
#' @noRd
plot.cov = function(mat, clamp=TRUE, palette=colryb,
                    breaks=NULL, breakscol='black', blty=1,zlim=NULL){
    if (clamp) mat = clamp.mat(mat)
    if (is.null(zlim)){
        image(mat, axes=F, col=palette, useRaster=TRUE)
    } else { 
        image(mat, axes=F, col=palette, zlim=zlim, useRaster=TRUE)
    }
    box(lwd=0.5)
    if (!is.null(breaks)){
        abline(v=breaks,lty=blty,lwd=.25, col=breakscol)
        abline(h=breaks,lty=blty,lwd=.25, col=breakscol)
    }
}

#' Adust color transparency for hex color
#'
#' @param x color
#' @param alpha transparency level
#' @return adjusted color
#' @noRd
tsp.col = function(x, alpha=0.5){
    rr = col2rgb(x)
    rgb(t(rr)/255, alpha=alpha)
}

#' Calculate the breaks from a series of clusters
#'
#' @param acut
#' @param start
#' @param end
#' @return breaks
#' @noRd
calc.breaks.acut <- function(acut, start=0, end=1){
    cuts <- cumsum(rle(acut)$lengths)
    step <- head(diff(seq(start,end,length.out=length(acut))),1)
    cuts <-  cuts[-length(cuts)] -1
    cuts * step + step/2
}

#' Calculate the breaks from the cluster object 
#'
#' @param ht cluster object
#' @param nclust number of clusters
#' @param cocl clusters
#' @return breaks
#' @noRd
calc.breaks <- function(ht, nclust, cocl){
    acut <- cutree(ht, nclust)[cocl]
    cuts <- cumsum(rle(acut)$lengths)
    step <- head(diff(seq(0,1,length.out=length(cocl))),1)
    cuts <-  cuts[-length(cuts)] -1
    cuts * step + step/2
}

#' Calculate the breaks and rectangles from two clusters
#'
#' @param hcls horizontal clusters
#' @param vcls vertical clusters
#' @param colset colors for the clusters
#' @param start
#' @param end
#' @return rectangles and colors
#' @noRd
calc.breaks.rect <- function(hcls, vcls, colset, start=0, end=1){
    # Calculate locations of breaks:
    vbreaks = calc.breaks.acut(vcls, start=start, end=end)
    hbreaks = calc.breaks.acut(hcls, start=start, end=end)
    # Matching cls:
    kv = unique(vcls)
    kh = unique(hcls)
    kb = sort(kh[kh %in% kv])
    vbreaks = c(par()$usr[1], vbreaks, par()$usr[2])
    hbreaks = c(par()$usr[3], hbreaks, par()$usr[4])
    names(hbreaks) = NULL
    names(vbreaks) = NULL
    # Match breaks to each other:
    rectdf = c()
    for (id in kb){
        ih = which(kh == id)
        iv = which(kv == id)
        rectdf = rbind(rectdf,
                       c(x1=vbreaks[iv], x2=vbreaks[iv + 1],
                         y1=hbreaks[ih], y2=hbreaks[ih + 1]))
    }
    rectdf = data.frame(rectdf)
    if (kb[1] == 0){
        vccols = c('grey', colset[kb[-1]])
    } else {
        vccols = colset[kb]
    }
    return(list(rectdf, vccols))
}

#' Diagonalize and order columns of centers of matrix:
#' 
#' NOTE: Old version that is incorrect, but maintained
#' @noRd
diag.mat = function(mat, ratio=0.5){
    ord = order(colSums(mat > 0.25)  > ratio * ncol(mat),
                apply(mat,2,which.max), decreasing=TRUE)
    mat = mat[,ord]
    cto = apply(mat,2,which.max)
    idx = colSums(mat > 0.25)  > ratio * ncol(mat)
    cto[idx] = 0
    return(list(mat, colnames(mat), cto))
}

#' Diagonalize and order columns of centers of matrix:
#' 
#' NOTE: Correct version
#' @noRd
diag.mat2 = function(mat, ratio=0.5, cutoff=0.25){
    ord = order(colSums(mat > cutoff)  > ratio * nrow(mat),
                apply(mat,2,which.max), decreasing=TRUE)
    mat = mat[,ord]
    cto = apply(mat,2,which.max)
    idx = colSums(mat > cutoff)  > ratio * nrow(mat)
    cto[idx] = 0
    return(list(mat, colnames(mat), cto))
}

#' Label runs in a column of labels
#' 
#' @param factor.labels Labels as factors
#' @param labels Labels
#' @param rdcol Colors for the labels
#' @noRd
label.runs <- function(factor.labels, labels, rdcol){
    rl = rle(as.numeric(factor.labels))
    cs = cumsum(rl$lengths)
    cb = c(0, cs[-length(cs)] + 1)
    loc = (cs + cb )/ 2 / length(labels)
    idx = which(rl$lengths > 2)
    lab = levels(labels)[rl$values[idx]]
    # Remove NONE and average consecutive values:
    kid = which(lab != 'NONE')
    col.lab = as.character(rdcol$COLOR[rl$values[idx[kid]]])
    lab = lab[kid]
    loci = loc[idx[kid]]
    rl = rle(lab)
    cs = cumsum(rl$lengths)
    cb = c(0, cs[-length(cs)] + 1)
    locfinal = c()
    lab = rl$values
    col.lab = col.lab[cs]
    locfinal = sapply(1:length(cs), function(i){mean(loci[cb[i]:cs[i]])})
    return(list(locfinal, lab, col.lab))
}

#' Process a distance matrix for the epimap samples 
#' 
#' @param mat
#' @noRd
process.dist.mat <- function(mat){
    if (dim(mat)[1] > 2){
        idx = grep("^BSS", colnames(mat))
        mat <- mat[idx, idx]
    }
    if (dim(mat)[1] > 2){
        rownames(mat) <- colnames(mat)
        # Fill in backwards:
        idx <- mat == 0
        mat[idx] <- t(mat)[idx]
        mat[is.na(mat)] <- 0
        mat[mat == 0] <- NA
        # Fill in missing with median guess: 
        medians1 = apply(mat, 1, mean, na.rm=T)
        medians2 = apply(mat, 2, mean, na.rm=T)
        repl = outer(medians1, medians2, '+') / 2
        mat[is.na(mat)] =  repl[is.na(mat)]
        diag(mat) = NA
        rownames(mat) <- colnames(mat)
        mat <- as.matrix(mat)
        mat.orig <- mat
        mat <- reord(mat)
        mat <- mat[,rownames(mat)]
    }
    return(mat)
}

#' Plot metadata image with appropriate colors
#' 
#' @param metamat Metadata matrix
#' @param colvals List of per-covariate colors for matrix
#' @noRd
meta.image <- function(metamat, colvals=NULL, labels=NULL, horiz=FALSE, cex=1, return.mat=FALSE, useRaster=FALSE){
    if (is.null(colvals)){
        colmat = metamat
    } else {
        colmat = metamat
        for (nam in colnames(colmat)){
            colmat[,nam] = colvals[[nam]][metamat[,nam]]
        }
    }
    if (is.null(labels)){ labels = sapply(colnames(colmat), capitalize) }
    cc = as.character(colmat)
    cc[is.na(cc)] <- 'white'
    cols = unique(cc)
    cf = factor(cc, levels=cols)
    nummat = matrix(as.numeric(cf), nrow=nrow(colmat))
    if (!return.mat){
        if (horiz==FALSE){
            image(nummat, axes=F, col=cols, useRaster=useRaster)
            if (cex > 0) {
                text(y=seq(0, 1, length.out=length(labels)),
                     x=par()$usr[1]-0.001*(par()$usr[2]-par()$usr[1]), 
                     labels=labels, srt=0, adj=1, xpd=TRUE,cex=cex)
            }
        } else {
            image(t(nummat), axes=F, col=cols, useRaster=useRaster)
            if (cex > 0){
                text(x=seq(0, 1, length.out=length(labels)),
                     y=par()$usr[3]-0.001*(par()$usr[4]-par()$usr[3]), 
                     labels=labels, srt=90, adj=1, xpd=TRUE,cex=cex)
            }
        }
    } else {
        return(list(nummat, cols))
    }
}

# ---------------
# Text functions:
# ---------------

#' Capitalize string
#' 
#' @param x string to capitalize
#' @noRd
capitalize <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1,1)), substring(s, 2),
      sep="", collapse=" ")
}

#' Split text by widths
#' 
#' @param x string to split
#' @param x text width
#' @noRd
split.text = function(x, width=90){
    x = as.character(x)
    sl = as.numeric(gregexpr(" ", x)[[1]])
    slat = c(1, sl[which(diff(sl %% width) < 0)], nchar(x))
    x2 = ""
    for (i in 2:length(slat)){
        x2 = paste0(x2, substr(x, slat[i-1], slat[i]), "\n")
    }
    x2 = gsub(" \n ", "\n", x2)
    x2 = gsub("\n$","",x2)
    return(x2)
}


#' Space labels out for 1d
#' 
#' @param xx positions of text
#' @param box.pad box size/padding for text
#' @param lim limits of space to use
#' @noRd
space.1d = function(xx, box.pad=0.02, lim=c(0,1)){
    # Identify centers of runs, shift things away from them:
    dx = round(diff(xx), 4)
    xl = dx < box.pad
    rl = rle(xl)
    rl$lengths
    rdf = data.frame(l=rl$l, v=rl$v, cs=cumsum(rl$l))
    rdf$cent = rdf$cs - (rdf$l - 1) / 2
    pivots = rdf$cent[rdf$v] + 1
    x2 = xx
    nx = length(xx)
    # Can't just do this, might perturb the order.
    for (pv in pivots){
        if (round(pv,0) == pv){
            bot = pv - 1
            top = pv + 1
        } else {
            mx = mean(dx[c(pv - .5, pv + .5)])
            dx[pv - .5] = mx - box.pad / 2
            if (pv + .5 < nx - 1){ dx[pv + .5] = mx + box.pad / 2}
            bot = pv - .5 - 1
            top = pv + .5 + 1
        }
        # Go through and fix ripples:
        for (i in rev(1:bot)){
            if (i >= 1){
                if (x2[i + 1] - x2[i] <= box.pad){
                    x2[i] = x2[i + 1] - box.pad
                } else { break }
            }
        } 
        for (i in top:(nx-1)){
            if (i > 1 & i <= nx){
                if (x2[i] - x2[i - 1] <= box.pad){
                    x2[i] = x2[i - 1] + box.pad
                } else { break }
            }
        }
    }
    x2[x2 > lim[2]] = lim[2]
    x2[x2 < lim[1]] = lim[1]
    # Fix after lim, from outside in:
    for (i in 1:round(nx/2)){
        if (x2[i + 1] - x2[i] <= box.pad){
            x2[i + 1] = x2[i] + box.pad
        } else { break }
    } 
    for (i in rev(round(nx/2):nx)){
        if (x2[i] - x2[i - 1] <= box.pad){
            x2[i - 1] = x2[i] - box.pad
        } else { break }
    } 
    return(x2)
}

#' V2 of Upset plot for set intersections 
#' 
#' requires ComplexHeatmap
#' @noRd
UpSet.v2 = function(m, comb_col = "black", set_col = 'black', 
                    lwd = 2, pt_size = unit(3, "mm"), 
                    bg_col = "#F0F0F0", bg_pt_col = "#CCCCCC",
                    set_order = order(set_size(m), decreasing = TRUE), 
                    comb_order = if (attr(m, "set_on_rows")) { 
                        order.comb_mat(m[set_order, ], decreasing = TRUE)
                    } else {
                        order.comb_mat(m[, set_order], decreasing = TRUE)
                    }, top_annotation = upset_top_annotation(m),
                    right_annotation = upset_right_annotation(m),
                    row_names_side = "left", ...){
    require(ComplexHeatmap)
    set_on_rows = attr(m, "set_on_rows")
    mode = attr(m, "mode")
    m2 = m
    class(m2) = "matrix"
    pt_size = pt_size
    lwd = lwd
    if (set_on_rows) {
        n_comb = ncol(m)
        n_sets = nrow(m)
        if (length(comb_col == 1)){ comb_col = rep(comb_col, n_comb) }
        if (length(set_col == 1)){ set_col = rep(set_col, n_sets) }
        layer_fun = function(j, i, x, y, w, h, fill) {
            nr = round(1/as.numeric(h[1]))
            nc = round(1/as.numeric(w[1]))
            subm = matrix(pindex(m2, i, j), nrow = nr, byrow = FALSE)
            for (k in seq_len(nr)) {
                if (k%%2) {
                    grid.rect(y = k/nr, height = 1/nr, just = "top", 
                              gp = gpar(fill=bg_col[1], col=NA))
                } else {
                    if (length(bg_col) > 1) {
                        grid.rect(y = k/nr, height = 1/nr, just = "top", 
                                  gp = gpar(fill=bg_col[2], col=NA))
                    }
                }
            }
            jj = unique(j)
            for (k in seq_len(nc)) {
                if (sum(subm[, k]) >= 2) {
                    i_min = min(which(subm[, k] > 0))
                    i_max = max(which(subm[, k] > 0))
                    grid.lines(c(k - 0.5, k - 0.5)/nc,
                               (nr - c(i_min, i_max) + 0.5)/nr,
                               gp = gpar(col = comb_col[jj[k]], lwd = lwd))
                }
            }
            grid.points(x, y, size = pt_size, pch = 16,
                        gp = gpar(col = ifelse(pindex(m2, i, j), 
                                               set_col[i], NA)))
        }
        # Bars on top:
        ra = top_annotation
        if (length(ra) == 1) {
            ta_call = substitute(top_annotation)
            ta_call = as.list(ta_call)
            if (as.character(ta_call[[1]]) == "upset_top_annotation") {
                if (!"gp" %in% names(as.list(ta_call))) {
                    ra@anno_list[[1]]@fun@var_env$gp$fill = comb_col
                    ra@anno_list[[1]]@fun@var_env$gp$col = comb_col
                }
            }
        }
        sa = right_annotation
        if (length(sa) == 1) {
            ta_call = substitute(right_annotation)
            ta_call = as.list(ta_call)
            if (as.character(ta_call[[1]]) == "upset_right_annotation") {
                if (!"gp" %in% names(as.list(ta_call))) {
                    sa@anno_list[[1]]@fun@var_env$gp$fill = set_col
                    sa@anno_list[[1]]@fun@var_env$gp$col = set_col
                }
            }
        }
        ht = Heatmap(m2, cluster_rows = FALSE, cluster_columns = FALSE, 
                     rect_gp = gpar(type = "none"), layer_fun = layer_fun, 
                     show_heatmap_legend = FALSE, top_annotation = ra, 
                     right_annotation = sa,
                     row_names_side = row_names_side, 
                     row_order = set_order, column_order = comb_order, 
                     ...)
    } else {
        n_comb = nrow(m)
        if (length(comb_col == 1)) 
            comb_col = rep(comb_col, n_comb)
        layer_fun = function(j, i, x, y, w, h, fill) {
            nr = round(1/as.numeric(h[1]))
            nc = round(1/as.numeric(w[1]))
            subm = matrix(pindex(m2, i, j), nrow = nr, byrow = FALSE)
            for (k in seq_len(nc)) {
                if (k %% 2) {
                    grid.rect(x = k/nc, width = 1/nc, just = "right", 
                              gp = gpar(fill=bg_col, col=NA))
                }
            }
            grid.points(x, y, size = pt_size, pch = 16,
                        gp = gpar(col = ifelse(pindex(m2, i, j),
                                               set_col[j], bg_pt_col)))
            # comb_col[i], "#CCCCCC")))
            ii = unique(i)
            for (k in seq_len(nr)) {
                if (sum(subm[k, ]) >= 2) {
                    i_min = min(which(subm[k, ] > 0))
                    i_max = max(which(subm[k, ] > 0))
                    grid.lines((c(i_min, i_max) - 0.5)/nc, (nr - c(k, k) + 0.5)/nr, gp = gpar(col = comb_col[ii[k]], lwd = lwd))
                }
            }
        }
        ra = right_annotation
        if (length(ra) == 1) {
            ta_call = substitute(top_annotation)
            ta_call = as.list(ta_call)
            if (as.character(ta_call[[1]]) == "upset_right_annotation") {
                if (!"gp" %in% names(as.list(ta_call))) {
                    ra@anno_list[[1]]@fun@var_env$gp$fill = comb_col
                    ra@anno_list[[1]]@fun@var_env$gp$col = comb_col
                }
            }
        }
        ht = Heatmap(m2, cluster_rows = FALSE, cluster_columns = FALSE, 
                     rect_gp = gpar(type = "none"), layer_fun = layer_fun, 
                     show_heatmap_legend = FALSE, 
                     top_annotation = top_annotation, right_annotation = ra,
                     row_order = comb_order, column_order = set_order, ...)
    }
    ht
}


#' Make a transformation matrix for factor data
#' 
#' @noRd
make.tform = function(x, norm=FALSE, u=NULL){
    if (is.null(u)){u = unique(x)}
    tmat = sapply(u, function(y){1 * (x == y)})
    if (norm){
        tmat = sweep(tmat, 2, apply(tmat, 2, sum),  '/')
    }
    return(tmat) }


#' Plot distance symmetric matrix
#' 
#' @param dt distance object
#' @noRd
plot.dt.sym = function(dt, nbreak=NULL, txtcol='black'){
    mat = as.matrix(dt)
    ht <- hclust(dt, method='ward.D')
    cocl <- order.optimal(dt, ht$merge)$order
    rn <- names(cocl)[cocl]
    mat = mat[rn,rn]
    colramp=rev(colred)
    colramp=colrb
    colramp=colspec
    diag(mat) = NA
    # Plotting function:
    par(mar=c(.5,10,.5,.5))
    par(yaxs='i')
    par(xaxs='i')
    minor = 1
    # image(mat, axes=F, zlim=c(0,1), col=colramp, useRaster=TRUE)
    image(mat, axes=F, col=colramp, useRaster=TRUE)
    yt = seq(0,1,length.out=ncol(mat))
    xt = par()$usr[1] - 0.01*diff(par()$usr[1:2])
    yaxlab=rownames(mat)
    if (length(txtcol) > 1){ txtcol = txtcol[cocl] }
    text(y=yt, x=xt, labels=yaxlab, srt=0, adj=1, xpd=TRUE,cex=1, col=txtcol)
    # Generate breaks:
    if (!(is.null(nbreak))){ 
        breaks = calc.breaks(ht, nbreak, cocl)
        abline(v=breaks,lty=2,lwd=.5)
        abline(h=breaks,lty=2,lwd=.5)
    }
}


#' Compute go enrichments by over-representation test
#' 
#' Requires clusterProfiler
#' @param dt distance object
#' @noRd
go.enr = function(x, ont='BP', tag=NULL, allx=gmdf$symbol,
                  minsize=10, maxsize=500){
    require(clusterProfiler)
    eg = bitr(allx, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
    eg$val = 1 * (eg$SYMBOL %in% x)
    genes = eg$ENTREZID[eg$val > 0]
    # Run enrichment
    gse.res = enrichGO(gene = genes, universe = eg$ENTREZID,
                       OrgDb = org.Hs.eg.db, ont = ont,
                       pAdjustMethod = "BH",
                       minGSSize=minsize, maxGSSize=maxsize,
                       pvalueCutoff = 0.1, readable = TRUE)
    gsedf = as.data.frame(gse.res)
    # Get top term for each combination:
    if (nrow(gsedf) > 0){
        u = unique(gsedf$geneID)
        df = c()
        for (set in u){
            sdf = gsedf[gsedf$geneID == set,]
            df = rbind(df, sdf[1,])
        }
        df = df[order(df$qvalue), ]
        df$Description = sapply(df$Description, width=90, split.text)
        df$ont = 'BP'
        if (!is.null(tag)){ df$tag = tag } else { df$tag = '' }
        out = df[, c('Description', 'qvalue', 'geneID', 'ont','tag')]
    } else {out = c() }
    return(out)
}
