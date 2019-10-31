#!/usr/bin/R
# ----------------------
# General use functions:
# ----------------------

#' Set domain and working directory to a specific project
#'
#' @export
set_proj = function(project){
    domain = system("hostname -d", intern=TRUE)
    if (length(domain) == 0) domain = ''
    if (domain == 'broadinstitute.org'){
        maindir=paste0('~/data/', project, '/')
    } else {
        maindir=paste0('~/', project, '/')
    }
    assign("dbdir", paste0(maindir, 'db/'), envir = .GlobalEnv)
    assign("bindir", paste0(maindir, 'bin/'), envir = .GlobalEnv)
    assign("img", paste0(maindir, 'img/'), envir = .GlobalEnv)
    options(scipen=45)
    setwd(dbdir)
}

#' List objects + their sizes:
#' 
#'  https://stackoverflow.com/questions/1358003/
#'
#' @export
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

#' Shorthand for list objects + their sizes:
#' 
#'  https://stackoverflow.com/questions/1358003/
#'
#' @export
lsos <- function(..., n=10) {
    .ls.objects(..., order.by="Size", decreasing=TRUE, head=TRUE, n=n)
}

#' Get percentiles for each value
#'
#' @param x
#' @return percentile rank
#' @export
perc.rank <- function(x){ trunc(rank(x))/length(x) }

#' Make a transformation matrix for factor data
#' 
#' @export
make.tform = function(x, norm=FALSE, u=NULL){
    if (is.null(u)){u = unique(x)}
    tmat = sapply(u, function(y){1 * (x == y)})
    if (norm){
        tmat = sweep(tmat, 2, apply(tmat, 2, sum),  '/')
    }
    return(tmat) 
}

#' Sigmoid function
#'
#' @export
sigm = function(z) {1 / (1 + exp(-z))}

#' Logit: Inverse of sigmoid function
#'
#' @export
rsigm = function(p) {log(p / (1 - p))}

#' Logit with scale factor
#'
#' @export
logit = function(x, scale=1){ x = x / scale; log(x / (1 - x)) }

