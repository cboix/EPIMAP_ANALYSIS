#!/usr/bin/R
# -------------------------------------
# Auxiliary functions for epimap paper:
# -------------------------------------

#' Label runs in a column of labels
#' 
#' @param factor.labels Labels as factors
#' @param labels Labels
#' @param rdcol Colors for the labels
#' @export
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
#' @export
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
#' @export
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

#' Upset plot (V2 for coloring bars) for set intersections 
#' 
#' Note: Requires ComplexHeatmap
#' @export
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

#' Plot distance symmetric matrix
#' 
#' @param dt distance object
#' @export
plot.dt.sym = function(dt, nbreak=NULL, txtcol='black'){

    require(cba)
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
#' @export
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
