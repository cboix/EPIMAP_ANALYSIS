#!/usr/bin/R
# General requisite libraries / functions: 
library(RColorBrewer)
library(reshape2)
library(plyr)
library(tidyr)
library(Matrix)
col <- colorRampPalette(brewer.pal(n=9,name="Blues"))(100) # For old plots
col1 <- colorRampPalette(brewer.pal(n=9,name="Blues"))(100)
col2 <- colorRampPalette(brewer.pal(n=7,name="Reds"))(100)
col3 <- colorRampPalette(brewer.pal(n=7,name="Purples"))(100)
col4 <- colorRampPalette(brewer.pal(n=7,name="Greens"))(100)
colrb <- colorRampPalette(brewer.pal(n=7,name="RdBu"))(100)
colrwb <- colorRampPalette(c('darkred','white','royalblue'))(200)
colramp <- colorRampPalette(c('white','royalblue'))(100)
colred <- colorRampPalette(c('white','darkred'))(100)
colspec <- colorRampPalette(brewer.pal(n=11,name="Spectral"))(100)
colryb <- colorRampPalette(brewer.pal(n=11,name="RdYlBu"))(100)
colryg <- colorRampPalette(brewer.pal(n=11,name="RdYlGn"))(100)
colpair <- colorRampPalette(brewer.pal(n=12,name="Paired"))(100)

options(scipen=45) # so we dont get issues writing integers into bedfiles

# Set directories: 
domain <- Sys.getenv('DOMAINNAME')
if (domain == 'broadinstitute.org'){ 
    setwd('/broad/compbio/cboix/EPIMAP_ANALYSIS/db/')
    img <- '/broad/compbio/cboix/EPIMAP_ANALYSIS/img/'
} else { 
    setwd('~/EPIMAP_ANALYSIS/db/')
    img <- '~/EPIMAP_ANALYSIS/img/'
    # Libraries not in broad:
    library(cba)
}

# Useful functions:
corner <- function(m,x=5){ m[1:x,1:x]}

# Reorder according to order.optimal (cba package)
reord <- function(mat, measure='eJaccard', method='complete') {
    # requires cba
    dt <- dist(mat, measure)
    ht <- hclust(dt, method=method)
    cocl <- order.optimal(dt, ht$merge)$order
    reord <- names(cocl)[cocl]
    mat[reord,]
}

get_optimal_ordering <- function(distmat, method='average') {
    library(cba)
    hc <- hclust(distmat, method=method)
    optim <- order.optimal(distmat, hc$merge)
    invisible(as.numeric(optim$order))
}


pivot.tomatrix = function(df, key, value){
    require(dplyr)
    wide = spread(df, key, value)
    mat = as.matrix(wide[,-1])
    rownames(mat) = wide[,1]
    return(mat)
}

# List objects + their sizes:
# https://stackoverflow.com/questions/1358003/
.ls.objects <- function (pos = 1, pattern, order.by,
                         decreasing=FALSE, head=FALSE, n=5) {
    napply <- function(names, fn) sapply(names, function(x)
                                         fn(get(x, pos = pos)))
    names <- ls(pos = pos, pattern = pattern)
    obj.class <- napply(names, function(x) as.character(class(x))[1])
    obj.mode <- napply(names, mode)
    obj.type <- ifelse(is.na(obj.class), obj.mode, obj.class)
    obj.size <- napply(names, object.size)
    obj.dim <- t(napply(names, function(x)
                        as.numeric(dim(x))[1:2]))
    vec <- is.na(obj.dim)[, 1] & (obj.type != "function")
    obj.dim[vec, 1] <- napply(names, length)[vec]
    out <- data.frame(obj.type, obj.size, obj.dim)
    names(out) <- c("Type", "Size", "Rows", "Columns")
    if (!missing(order.by))
        out <- out[order(out[[order.by]], decreasing=decreasing), ]
    if (head)
        out <- head(out, n)
    out
}
# shorthand
lsos <- function(..., n=10) {
    .ls.objects(..., order.by="Size", decreasing=TRUE, head=TRUE, n=n)
}
