#!/usr/bin/env Rscript
# -----------------------------
# Plot parsed raw GREAT outputs
# Adapted from WM
# -----------------------------
domain = system("hostname -d", intern=TRUE)
if (domain == 'broadinstitute.org'){
    bindir='~/data/EPIMAP_ANALYSIS/bin/'
} else {
    bindir='~/EPIMAP_ANALYSIS/bin/'
}
# Load functions:
# library(colormap)
source(paste0(bindir,'general_EPIMAP_ANALYSIS.R'))
source(paste0(bindir, 'great_enrichment/general_WM_functions.R'))
source(paste0(bindir, 'great_enrichment/auxfunctions_GREAT.R'))

args=(commandArgs(TRUE))
if (length(args)==0) {
    stop("No arguments supplied.")
} else {
    infofile = args[1]
    resultsdir = args[2]
    filetype = args[3]
}

# resultsdir='/home/cboix/go_enrichments/'
# resultsdir='./'
# filetype='png'

# Directories:
parseddir <- paste(resultsdir, "parsed_results", sep="/");
outputdir <- paste(resultsdir, "figures_aggregate", sep="/");
dir.create(outputdir, showWarnings=FALSE, recursive=TRUE);

# info = read.delim(infofile, header=F) 
# names(info) = c('file', 'name')

# source("../../code/general.R")
# #source("../general.R")
# source("general.R");
# source("../../Roadmap_25states/general.R")

noreorder = FALSE  # TRUE

# Set numeric plotting arguments :
# Heuristics...
n <- 10;
cexRow <- (1/n)*40;
marginRow <- (1/sqrt(n))*120
# k <- length(epitopes);
k = nrow(info)

# Load data:
load(paste(parseddir, "/res_all_per_onto.RData", sep=""))
ontos <- names(res_all)  # Ontology names

# For each ontology, plot a heatmap (no abbreviating/raw):
for (onto in ontos) {
    print(onto)
    res <- res_all[[onto]];
    resmat <- as.matrix(res, ncol=k);
    resmat = reorder_mat(res, noreorder=noreorder)
    nam <- gsub(" ", "_", onto);
    if (class(resmat) == "matrix") {
        plot_heatmap_image(resmat, paste(nam, sep="_"), rowScale=12, filetype=filetype, dims=c(8,11));
        plot_heatmap_horiz_image(resmat, paste(nam, sep="_"), filetype=filetype, dims=c(11,8));
    }
}


# Abbreviate the names of the rows:
abbrev <- sapply(names(res_all), function(x) {
  paste(sapply(strsplit(x, ""), function(y) { 
    sel <- which(y %in% LETTERS);
    sel <- sort(unique(c(sel, sel+1)));
    sel <- sel[y[sel] %in% LETTERS | y[sel] %in% letters]
    y[sel]
  }), collapse="")
})

rownames_filt <- as.character(unlist(sapply(names(res_all), function(x) { paste(abbrev[x], rownames(res_all[[x]]), sep=" - ") })));

# Join all ontologies into a complete matrix:
res_all_concat <- do.call(rbind, res_all)
resmat <- as.matrix(res_all_concat, ncol=k);

# Order the complete matrix by binary 
# distance of significant results.
resmat_row_ord <- get_optimal_ordering(dist(resmat > 2,
                                            method="binary"));
if (noreorder){
    resmat_col_ord <- 1:ncol(resmat)
} else {
    resmat_col_ord <- get_optimal_ordering(dist(t(resmat > 2),
                                                method="binary"));
}
resmat <- resmat[resmat_row_ord, resmat_col_ord];

# Save reordered + selected rownames:
rownames_filt <- rownames_filt[resmat_row_ord];
save(resmat, rownames_filt, file="resmat_rownames_filt.RData")

NWORDS=50
wordlets <- get_summary_terms(rownames(resmat), NWORDS)
wordlets_pos <- seq(1, nrow(resmat), length.out=length(wordlets)+1);
wordlets_pos <- wordlets_pos[-1] - diff(wordlets_pos)/2;

# Get dist matrix of reordered matrix + save:
distmat <- dist(t(resmat > 2), method="binary");
save(distmat, file="distmat_GREAT.RData")

# distmat <- dist(t(resmat))
hc <- hclust(distmat, method="average");
optim <- order.optimal(distmat, hc$merge);
hc$order <- optim$order
hc$merge <- optim$merge;
  
cols=c("white", "yellow", "orange", "red");
cutres <- cutree(hc, k=21)
cutres[cutres %in% which(tabulate(cutres) == 1)] <- 0;
newidx <- as.numeric(factor(as.character(cutres)));
# rowcols <- c("white", brewer.pal(max(newidx)-1, "Set1"))[newidx];
setcolors = c(brewer.pal(9,'Set1'),brewer.pal(8,'Set2'),brewer.pal(12,'Set3'))
rowcols <- c("white", setcolors)[newidx];
  
plotwidth=nrow(resmat)/400;
cexRow <- 1/log10(ncol(resmat)) * 1.5;


# Heatmaps with dendrograms:
for (i in 1:length(abbrev)) {
    print(names(abbrev)[i])
    wordlets_sub <- get_summary_terms(sapply(strsplit(rownames_filt, " - "), function(y) if (y[1] == abbrev[i]) y[2] else NA), NWORDS);
    plotfile(paste(outputdir, "/heatmap_pvals_horiz_all_w_dendro_annot_", abbrev[i], sep=""), height=7, width=plotwidth, type=filetype)
    heatmap(t(resmat), Rowv=as.dendrogram(hc), Colv=NA, col=cols,
            breaks=c(0, 2, 3, 4, max(c(5, max(resmat[is.finite(resmat)])+1))),
            main=paste("Annotated using", names(abbrev)[i]),
            scale="none",
            margins=c(24, 8),
            labCol=NA,
            cexRow=cexRow,
            RowSideColors=rowcols,
            cex.main=0.8,
            add.expr={box(lwd=1); axis(side=1, at=wordlets_pos,
                                       labels=wordlets_sub, las=2,
                                       cex.axis=1/log10(length(wordlets))*1.5,
                                       tick=FALSE, line=-0.5)});
    dev.off()
} 


# Heatmaps without dendrograms:
for (i in 1:length(abbrev)) {
    print(names(abbrev)[i])
    wordlets_sub <- get_summary_terms(sapply(strsplit(rownames_filt, " - "), function(y) if (y[1] == abbrev[i]) y[2] else NA), 40);
    plotfile(paste(outputdir, "/heatmap_pvals_horiz_all_w_annot_", abbrev[i], sep=""), height=7, width=plotwidth, type=filetype)
    heatmap(t(resmat), Rowv=NA, Colv=NA, col=cols,
            breaks=c(0, 2, 3, 4, max(c(5, max(resmat[is.finite(resmat)])+1))),
            main=paste("Annotated using", names(abbrev)[i]),
            scale="none",
            margins=c(24, 8),
            labCol=NA,
            cexRow=cexRow,
            RowSideColors=rowcols,
            cex.main=0.8,
            add.expr={box(lwd=1); axis(side=1, at=wordlets_pos,
                                       labels=wordlets_sub, las=2,
                                       cex.axis=1/log10(length(wordlets))*1.5,
                                       tick=FALSE, line=-0.5)});
    dev.off()
} 




# Heatmaps without dendrograms:
print(names(abbrev)[i])
plotfile(paste(outputdir, "/heatmap_pvals_horiz_all_w_annot_full", sep=""), height=7, width=plotwidth, type=filetype)

heatmap(t(resmat), Rowv=NA, Colv=NA, col=cols,
        breaks=c(0, 2, 3, 4, max(c(5, max(resmat[is.finite(resmat)])+1))),
        main=paste("Annotated using All"),
        scale="none",
        margins=c(24, 8),
        labCol=NA,
        cexRow=cexRow,
        RowSideColors=rowcols,
        cex.main=0.8,
        add.expr={box(lwd=1); axis(side=1, at=wordlets_pos,
                                   labels=wordlets, las=2,
                                   cex.axis=1/log10(length(wordlets))*1.5,
                                   tick=FALSE, line=-0.5)});

dev.off()

