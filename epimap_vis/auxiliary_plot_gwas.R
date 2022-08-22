#!/usr/bin/R
# ------------------------------------
# Functions for plotting GWAS results:
# ------------------------------------
library(RColorBrewer)
library(cba)
# library(scales)
# library(jsonlite)
# library(reshape2)
# library(httr)
# library(ape)
# library(infotheo)
# library(dendextend)

# Colors:
spec <- colorRampPalette(brewer.pal(n=9,name="Spectral"))(100)
rdgy <- colorRampPalette(brewer.pal(n=11,name="RdGy"))(100)
colr <- colorRampPalette(brewer.pal(n=9,name="Reds"))(100)
colb <- colorRampPalette(brewer.pal(n=9,name="Blues"))(100)


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


# Reorder according to order.optimal (cba package)
reord <- function(mat, measure='eJaccard', method='complete') {
    # requires cba
    dt <- dist(mat, measure)
    ht <- hclust(dt, method=method)
    cocl <- order.optimal(dt, ht$merge)$order
    reord <- names(cocl)[cocl]
    mat[reord,]
}

#' Diagonalize and order columns of centers of matrix:
#' 
#' NOTE: Correct version
#' @noRd
diag.mat = function(mat, ratio=0.5, cutoff=0.25){
    ord = order(colSums(mat > cutoff)  > ratio * nrow(mat),
                apply(mat,2,which.max), decreasing=TRUE)
    mat = mat[,ord]
    cto = apply(mat,2,which.max)
    idx = colSums(mat > cutoff)  > ratio * nrow(mat)
    cto[idx] = 0
    return(list(mat, colnames(mat), cto))
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


parpos = function(axis, shift){
    # NOTE: par()$usr is x1, x2, y1, y2
    if (axis == 1) { # X axis 
        par()$usr[1] - shift * diff(par()$usr[1:2])
    } else { # Y axis
        par()$usr[3] - shift * diff(par()$usr[3:4])
    }
}


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

process.gwas.matrix = function(mat, meta, odf, rdcol){
    # Cluster first:
    sampord = rownames(reord(t(mat) > 3, 'Jaccard'))
    studyord = rownames(reord(mat > 3, 'eJaccard'))
    mat = mat[studyord, sampord]

    # Order first by group, then according to odf:
    rnmat = colnames(mat)
    labelsmat = meta[rnmat, 'GROUP']
    labelsmat = as.character(labelsmat)
    labelsmat = factor(labelsmat, levels=odf$GROUP)
    rngroup = rnmat[order(labelsmat, decreasing=TRUE)]

    # Params for plotting:
    labels = meta[rngroup, 'GROUP']
    faclabels = as.matrix(as.numeric(labels))
    acutgroup = as.numeric(labels)
    dist.breaks = calc.breaks.acut(acutgroup)
    names(acutgroup) = rngroup
    lablist = label.runs(faclabels, labels, rdcol)
    colset = as.character(rdcol$COLOR)

    # Diagonalize, ignore all <2
    tmat = t(mat)
    ll = diag.mat(tmat[rngroup,], ratio=.25)
    tmp = ll[[1]]
    cto = ll[[3]] # For breaks
    diagord = colnames(tmp)
    vcut =  c(cto[cto ==0], acutgroup[cto])
    vbreaks = calc.breaks.acut(vcut)
    return(list(rngroup = rngroup,
                diagord = diagord,
                labels = labels,
                acut = acutgroup,
                vcut = vcut,
                colset = colset,
                lablist = lablist))
}


plot.gwas.matrix.flipped <- function(mat, attr, colvals, meta, metamat, rdcol,
                             highlight=NULL, palette=colb, zmax=20, 
                             subset.samples=NULL, subset.gwas=NULL){
    # Subset matrix if necessary:
    diagord = attr$diagord
    rngroup = attr$rngroup
    lablist = attr$lablist
    acut = attr$acut
    vcut = attr$vcut
    if (!is.null(subset.samples)){
        rind = rngroup %in% subset.samples
        rngroup = rngroup[rind] 
        acut = acut[rind]
        hbreaks = NULL
        # TODO: should further restrict
        labels = meta[rngroup, 'GROUP']
        faclabels = as.matrix(as.numeric(labels))
        print(length(labels))
        lablist = label.runs(faclabels, labels, rdcol)
        print(lablist)
    }
    hbreaks = calc.breaks.acut(acut)
    if (!is.null(subset.gwas)){
        gind = diagord %in% subset.gwas
        diagord = diagord[gind]
        vcut = vcut[gind]
        gwascex=0.85
    } else {gwascex = 0.5}
    vbreaks = calc.breaks.acut(vcut) 
    # TODO: Throw error if nothing in matrix?
    mat[mat > zmax] = zmax
    mat = t(mat[rev(diagord), rngroup])
    # TODO: add subsets here:
    layout(matrix(c(1,2),2,1), heights=c(1.2,10), TRUE)
    par(yaxs="i")
    par(xaxs="i")
    par(xpd=FALSE)
    lsp = 25
    par(mar=c(0, lsp, 12, 0.25))
    meta.image(metamat[rngroup,], colvals=colvals, cex=0, horiz=F, useRaster=TRUE)
    abline(h=par()$usr[3:4],lty=1,lw=0.5)
    abline(v=par()$usr[1:2],lty=1,lw=0.5)
    metaclass = sapply(colnames(metamat), capitalize)
    text(y=seq(0,1, length.out=ncol(metamat)),
         x=par()$usr[1]-0.001*(par()$usr[2]-par()$usr[1]), 
         labels=metaclass, srt=0, adj=1, xpd=TRUE, cex=.6)
    # Plot labels:
    box.pad = 0.02
    xw = lablist[[1]]
    if (length(xw) > 1){
        xx = space.1d(xw, box.pad=box.pad)
        xx = space.1d(xx, box.pad=box.pad)
    } else {xx = xw}
    rx = c(0.05, 0.15, 0.4, 0.45)
    x = par()$usr[4]+rx*(diff(par()$usr[3:4]))
    text(x=xx, y=x[4], labels=lablist[[2]], col=lablist[[3]],
         srt=90, adj=0, xpd=TRUE, cex=.8)
    par(xpd=TRUE)
    segments(x0=xw, y0=x[2], x1=xx, y1=x[3], col=lablist[[3]])
    segments(x0=xw, y0=x[1], x1=xw, y1=x[2], col=lablist[[3]])
    par(xpd=FALSE)
    # Legend on top:
    metanams = unlist(sapply(colnames(metamat)[1:4], function(x){paste0(capitalize(x),": ", names(colvals[[x]]))}))
    metacols = unlist(sapply(colnames(metamat)[1:4], function(x){colvals[[x]]}))
    legend("top", inset=c(0,-3.75), title="", metanams, col=metacols, 
           xpd=TRUE, ncol=4, pch=15, pt.cex=1.25, bty='n')

    par(mar=c(.25, lsp, .25, 0.25))
    image(mat, axes=FALSE,col=palette, zlim=c(0, zmax), useRaster=TRUE)
    abline(h=par()$usr[3:4],lty=1,lw=0.5)
    abline(v=par()$usr[1:2],lty=1,lw=0.5)
    if (!is.null(highlight)){
        txtcol = ifelse(rev(diagord) %in% highlight, 'red','black')
    } else { txtcol = 'black' }
    text(y=seq(0,1, length.out=length(diagord)), x=parpos(1, 0.005), 
         labels=rev(diagord), srt=0, adj=1, xpd=TRUE,cex=gwascex, col=txtcol)
    # Add rectangles
    rll = calc.breaks.rect(hcls=rev(vcut), vcls=acut, attr$colset)
    rectdf = rll[[1]]
    vccols = rll[[2]]
    rord = order(rectdf$y1)
    rectdf = rectdf[rord,]
    vccols = vccols[rord]
    abline(h=unique(c(rectdf$y1, rectdf$y2)),lty='dotted', 
           lwd=.5, col='darkgrey')
    abline(v=hbreaks,lty='dotted', lwd=.5, col='darkgrey')
    # Add grey - ubq rectangle:
    rectdf = rbind(c(x1=par()$usr[1], x2=par()$usr[2],
                     y1=rectdf$y2[nrow(rectdf)], y2=par()$usr[4]), rectdf)
    vccols = c('grey',vccols)
    rect(xleft=rectdf$x1, xright=rectdf$x2,
         ybottom=rectdf$y1, ytop=rectdf$y2,
         border=vccols, lwd=1)
    par(xpd=NA)
    rect(ybottom=rectdf$y1, ytop=rectdf$y2, xright=par()$usr[1], 
         xleft=par()$usr[1]-0.004*(par()$usr[2]-par()$usr[1]), 
         border='white', col=vccols, lwd=.25)
}


