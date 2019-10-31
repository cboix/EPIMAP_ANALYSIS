#!/usr/bin/env Rscript
# -----------------------------
# Aux functions for GO plots
# All functions from WM
# -----------------------------
library(caTools)
library(gtools)
library(GO.db)
library(RColorBrewer)
# library(gplots)
library(cba)
options(scipen=20)
options(warn=2)


# Define functions:
plot_heatmap_image <- function(resmat, nam, rowScale=12, 
                               filetype="pdf", dims=NULL,
                               cols=c("white", "yellow", "orange", "red")) {
    # Heuristics...
    nr <- nrow(resmat);
    nc <- ncol(resmat);
    if (is.null(dims)){
        height = 11
        width = height * nc / nr / rowScale
    } else {
        width = dims[1]
        height = dims[2]
    }
    if (!is.null(nr)) {
        cexRow <- 1/log10(nr) * 0.8
        cexCol <- 1/log10(nc) * 0.8
        resmat[!is.finite(resmat)] <- max(resmat[is.finite(resmat)])+10;
        plotfile(paste(outputdir, "/heatmap_pvals_vert_", nam, sep=""), height=height, width=width, type=filetype)
        par(mar=c(30,.25,2.5,.25))
        image(resmat, col=cols, axes=F)
        text(x=seq(0,1, length.out=nrow(resmat)),
             y=par()$usr[3]-0.01*(par()$usr[4]-par()$usr[3]), 
             labels=rownames(resmat), srt=90, adj=1, xpd=TRUE,cex=cexRow)
        box()
        mtext(gsub("_", " ", nam), side=3, line=.5, cex=cexRow * 2)
        dev.off()
    }
}

# Define functions:
plot_heatmap_horiz_image <- function(resmat, nam, rowScale=12, 
                                     filetype="pdf", dims=NULL,
                                     cols=c("white", "yellow", "orange", "red")) {
    # Heuristics...
    resmat = t(resmat)
    nr <- nrow(resmat);
    nc <- ncol(resmat);
    if (is.null(dims)){
        width = 11
        height = width * nr / nc / rowScale
    } else {
        width = dims[1]
        height = dims[2]
    }
    if (!is.null(nr)) {
        cexRow <- 1/log10(nr)
        cexCol <- 1/log10(nc) 
        resmat[!is.finite(resmat)] <- max(resmat[is.finite(resmat)])+10;
        plotfile(paste(outputdir, "/heatmap_pvals_horiz_", nam, sep=""), height=height, width=width, type=filetype)
        par(mar=c(.25,30,2.5,.25))
        image(resmat, col=cols, axes=F)
        text(y=seq(0,1, length.out=ncol(resmat)),
             x=par()$usr[1]-0.01*(par()$usr[2]-par()$usr[1]), 
             labels=colnames(resmat), srt=0, adj=1, xpd=TRUE,cex=cexRow)
        box()
        mtext(gsub("_", " ", nam), side=3, line=.5, cex=cexRow * 2)
        dev.off()
    }
}


# Define functions:
plot_heatmap <- function(resmat, nam, rowScale=12, filetype="pdf", 
                         cols=c("white", "yellow", "orange", "red")) {
    # Heuristics...
    nr <- nrow(resmat);
    nc <- ncol(resmat);
    if (!is.null(nr)) {
        cexRow <- 1/log10(nr) * 0.8
        cexCol <- 1/log10(nc) * 0.8
        resmat[!is.finite(resmat)] <- max(resmat[is.finite(resmat)])+10;
        plotfile(paste(outputdir, "/heatmap_pvals_vert_", nam, sep=""), height=nr/rowScale, width=7, type=filetype)
        try(heatmap.scale(resmat, width=8, height=nr/rowScale, Rowv=NA, Colv=NA, 
                          col=cols, breaks=c(0,2,3,4,max(c(5, max(resmat[is.finite(resmat)])+1))), 
                          scale="none", labRow=rownames(resmat), margins=c(10,12),
                          cexRow=cexRow, cexCol=cexCol, cex.main=cexCol, main=gsub("_", " ", nam),
                          add.expr={box(lwd=1);}));
        dev.off()
    }
}


plot_heatmap_horiz <- function(resmat, nam, filetype="pdf", cols=c("white", "yellow", "orange", "red")) {  
    nr <- nrow(resmat);
    nc <- ncol(resmat);
    cexRow <- 1/log10(nr)
    cexCol <- 1/log10(nc)
    plotfile(paste(outputdir, "/heatmap_pvals_horiz_", nam, sep=""), height=7, width=10, type="pdf")
    heatmap.scale(t(resmat), height=7, width=10, Rowv=NA, Colv=NA, 
                  col=cols, breaks=c(0,2,3,4,max(c(5, max(resmat[is.finite(resmat)])+1))), 
                  scale="none", margins=c(24,6), cexRow=cexCol, cexCol=cexRow, cex.main=0.8, main=gsub("_", " ", nam),
                  add.expr={box(lwd=1);});
    dev.off()
}

get_summary_terms <- function(dat, binsize=50, filtering=TRUE, tab_all=NULL) {
    terms_sel <- c();
    term_words <- strsplit(dat, "[ _,.]");
    if (is.null(tab_all)){
        tab_all <- sort(table(unlist(term_words)));
    }
    nbins <- floor(length(dat)/binsize);
    for (bin in 1:nbins) {
        i <- ((bin-1)*binsize)+1;
        j <- bin*binsize;
        tab <- table(unlist(term_words[i:j]));
        ratio <- log2((tab/binsize) / (tab_all[names(tab)]/length(dat)));
        ranks <- rank(ratio);
        if (filtering) ranks[tab < 2] <- NA;
        term_score <- sapply(term_words[i:j], function(x) {
                                 mean(ranks[names(ranks) %in% x], na.rm=T);
    });
        val <- which.max(term_score);
        terms_sel <- c(terms_sel, ifelse(length(val)==0, "", dat[i:j][val]));
    }
    invisible(terms_sel);
}

reorder_mat <- function(res, noreorder=FALSE, method='binary'){
    resmat <- as.matrix(res, ncol=k);
    hitnames <- rownames(res);
    if (length(hitnames) > 2) {
        if(noreorder){
            resmat_ord <- 1:ncol(resmat)
        } else {
            # TODO: May be better as eJaccard:
            resmat_ord <- get_optimal_ordering(dist(t(resmat > 2), method="binary"));
        }
        resmat <- resmat[,resmat_ord];
    }
    if (nrow(resmat) > 2) {
        resmat_ord <- get_optimal_ordering(dist(resmat > 2, method="binary"));
    } else {
        resmat_ord <- 1:nrow(resmat);
    }
    resmat <- resmat[resmat_ord,];
}
