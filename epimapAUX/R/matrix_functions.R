#!/usr/bin/R
# ---------------------------------------
# Matrix manipulation/analysis functions:
# ---------------------------------------

#' Show the top corner of a matrix:
#'
#' @param m matrix
#' @param x number of top components to show
#' @export
corner <- function(m,x=5){ m[1:x,1:x]}

#' Reorder matrix according to order.optimal (cba)
#'
#' @param mat matrix to reorder
#' @param measure distance measure from dist
#' @param method clustering method from hclust
#' @export
reord <- function(mat, measure='eJaccard', method='complete') {
    require(cba)
    dt <- dist(mat, measure)
    ht <- hclust(dt, method=method)
    cocl <- order.optimal(dt, ht$merge)$order
    reord <- names(cocl)[cocl]
    mat[reord,]
}

#' Reorder dist. matrix according to order.optimal (cba)
#'
#' @param distmat distance matrix
#' @param method clustering method from hclust
#' @export
get_optimal_ordering <- function(distmat, method='average') {
    require(cba)
    hc <- hclust(distmat, method=method)
    optim <- order.optimal(distmat, hc$merge)
    invisible(as.numeric(optim$order))
}

#' Turn a dataframe with key and value into a wide matrix
#'
#' NOTE: From 3-column df, where column 1 is another axis of interest
#'
#' @param df 3-column df
#' @param key column for pivot value
#' @param value column for values to put in matrix
#' @export
pivot.tomatrix = function(df, key, value){
    require(dplyr)
    wide = spread(df, key, value)
    mat = as.matrix(wide[,-1])
    rownames(mat) = wide[,1]
    return(mat)
}

#' Clamp a matrix at certain quantiles
#'
#' @param mat
#' @param alpha optional cutoff.
#' @return clamped matrix
#' @export
clamp.mat <- function(mat, alpha=0.01){
    quant = quantile(mat, c(0.01, 1 - 0.01), na.rm=T)
    mat[mat < quant[1]] = quant[1]
    mat[mat > quant[2]] = quant[2]
    return(mat)
}

#' Adust color transparency for hex color
#'
#' @param x color
#' @param alpha transparency level
#' @return adjusted color
#' @export
tsp.col = function(x, alpha=0.5){
    rr = col2rgb(x)
    rgb(t(rr)/255, alpha=alpha)
}

#' Plot a covariance matrix
#'
#' @param mat symmetric matrix
#' @param clamp whether to clamp
#' @export
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

#' Calculate the breaks from a series of clusters
#'
#' @param acut
#' @param start
#' @param end
#' @return breaks
#' @export
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
#' @export
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
#' @export
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
#' @export
diag.mat = function(mat, ratio=0.5){
    # Here, should have been nrow(mat):
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
#' @export
diag.mat2 = function(mat, ratio=0.5, cutoff=0.25){
    ord = order(colSums(mat > cutoff)  > ratio * nrow(mat),
                apply(mat,2,which.max), decreasing=TRUE)
    mat = mat[,ord]
    cto = apply(mat,2,which.max)
    idx = colSums(mat > cutoff)  > ratio * nrow(mat)
    cto[idx] = 0
    return(list(mat, colnames(mat), cto))
}


#' Faster jaccard distance
#' 
#' Calculate jaccard distance (1 - nint / nunion)
#' 
#' @param x First matrix for jaccard distance
#' @param y Second matrix. If NULL, uses transpose of x
#' @export
jacc.dist = function(x, y=NULL, verbose=TRUE){
    require(Matrix) # For if sparse
    if (is.null(y)){ 
        if (verbose){ print("Using t(x) as second matrix") }
        y = t(x) 
    }
    rsx = apply(x, 1, sum)
    rsy = apply(y, 2, sum)
    NTX = length(rsx)
    NTY = length(rsy)
    rsx = matrix(rep(rsx, NTX), nrow=NTX,ncol=NTX, byrow=T)
    rsy = matrix(rep(rsy, NTY), nrow=NTY,ncol=NTY, byrow=T)
    rsm = t(rsx) + rsy
    nint = x %*% y
    nunion = rsm - nint
    dist = 1.0 - nint / nunion 
    return(dist)
}


#' Cosine dist on two matrices:
#' 
#' @param x First matrix for jaccard distance
#' @param y Second matrix. If NULL, uses transpose of x
#' @export
cosine.dist = function(x, y=NULL, verbose=TRUE){
    if (is.null(y)){ 
        if (verbose){ print("Using t(x) as second matrix") }
        y = t(x) 
    }
    denom = sqrt(rowSums(x^2) %*% t(colSums(y^2)))
    dt = 1 - (x %*% y) / denom
    return(dt)
}
