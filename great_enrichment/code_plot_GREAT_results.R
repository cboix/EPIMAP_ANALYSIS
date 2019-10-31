source("../../code/general.R")
#source("../general.R")
source("general.R");
source("../../Roadmap_25states/general.R")
library(caTools)
library(gtools)
library(GO.db);
library(RColorBrewer)
library(gplots)

##########################################################################################################################

args=(commandArgs(TRUE))
if (length(args)==0) {
  stop("No arguments supplied.")
} else {
  eval(parse(text=args[[1]])) # parse first argument: cell_type
  eval(parse(text=args[[2]])) # parse first argument: filetype
}

resmat_hardcoded <- paste("f", c(83, 85, 99, 34, 24, 35, 18, 33, 81, 66, 71, 86, 41, 63, 27, 77, 61, 4, 10, 29, 64, 21, 74, 26, 57, 92, 
  97, 38, 28, 52, 87, 37, 15, 25, 51, 44, 98, 16, 75, 89, 17, 46, 22, 78, 20, 32, 65, 14, 40, 19, 42, 
  11, 50, 55, 82, 48, 36, 69, 6, 100, 31, 56, 54, 90, 5, 70, 73, 96, 12, 45, 72, 84, 60, 67, 47, 53, 
  62, 13, 79, 80, 76, 49, 59, 88, 91, 94, 23, 39, 9, 3, 2, 7, 43, 8, 30, 1, 93, 58, 68, 95), sep="");


##########################################################################################################################
### Define functions.
plot_heatmap <- function(resmat, nam, rowScale=12, filetype="pdf", cols=c("white", "yellow", "orange", "red")) {  
  # Heuristics...
  nr <- nrow(resmat);
  nc <- ncol(resmat);
  if (!is.null(nr)) {
    #cexRow <- 1/log10(nr) * 0.3
    cexRow <- 1/log10(nr) * 0.8
    cexCol <- 1/log10(nc) * 0.8
    resmat[!is.finite(resmat)] <- max(resmat[is.finite(resmat)])+10;
    plotfile(paste(outputdir, "/heatmap_pvals_vert_", nam, sep=""), height=nr/rowScale, width=7, type=filetype)
    try(heatmap.scale(resmat, width=8, height=nr/rowScale, Rowv=NA, Colv=NA, 
          col=cols, breaks=c(0,2,3,4,max(c(5, max(resmat[is.finite(resmat)])+1))), 
          scale="none", labRow=rownames(resmat), margins=c(10,12),
	      cexRow=cexRow, cexCol=cexCol, cex.main=cexCol, main=gsub("_", " ", nam),
          add.expr={box(lwd=1);}));
#          add.expr={box(lwd=1); abline(v=cumsum(bnds)+0.5, lwd=1)}));
#          add.expr={grid(ncol(resmat), nrow(resmat), lty=1, col=1, lwd=0.1); box(lwd=0.1);}));
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
#    add.expr={box(lwd=1); abline(v=cumsum(bnds)+0.5, lwd=1)}));
#    add.expr={grid(ncol(resmat), nrow(resmat), lty=1, col=1, lwd=0.1); box(lwd=0.1);}));
  dev.off()
}

get_summary_terms <- function(dat, binsize=50, filtering=TRUE) {
  terms_sel <- c();
  term_words <- strsplit(dat, "[ _,.]");
  tab_all <- sort(table(unlist(term_words)));

  nbins <- floor(length(dat)/binsize);
  for (bin in 1:nbins) {
    i <- ((bin-1)*binsize)+1;
    j <- bin*binsize;

    tab <- table(unlist(term_words[i:j]));
    ratio <- log2((tab/binsize) / (tab_all[names(tab)]/length(dat)));
    #ratio[tab < 3] <- NA;
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

##########################################################################################################################
  
# Set variables.
resultsdir <- paste("parsed_results", cell_type, sep="/");
outputdir <- paste("figures_aggregate", cell_type, sep="/");
dir.create(outputdir, showWarnings=FALSE, recursive=TRUE);

#########################################################################################################################

# Heuristics...
n <- 10; # dunno
cexRow <- (1/n)*40;
marginRow <- (1/sqrt(n))*120

#comparisons <- read.table("../WM20141020_comparative_epigenomics/comparisons.txt",
k <- length(epitopes);

#########################################################################################################################


load(paste(resultsdir, "/res_all_per_onto.RData", sep=""))
ontos <- names(res_all)

for (onto in ontos) {
  res <- res_all[[onto]];
  hitnames <- rownames(res);
  resmat <- as.matrix(res, ncol=k);

  if (length(hitnames) > 2) {
##    res_clust <- hclust(dist(t(resmat)));
##    resmat <- resmat[,res_clust$order];
    if(exists("resmat_hardcoded", envir=.GlobalEnv)) {
      resmat_ord <- resmat_hardcoded 
    } else {
      resmat_ord <- get_optimal_ordering(dist(t(resmat > 2), method="binary"));
    }
    resmat <- resmat[,resmat_ord];
  }

  if (nrow(resmat) > 2) {
    resmat_ord <- get_optimal_ordering(dist(resmat > 2, method="binary"));
    #resmat_ord <- order(apply(resmat, 1, which.max), decreasing=TRUE)
  } else {
    resmat_ord <- 1:nrow(resmat);
  }
  resmat <- resmat[resmat_ord,];

  nam <- gsub(" ", "_", onto);
#  plot_heatmap(resmat, nam, filetype=filetype);
#  plot_heatmap(resmat, paste(nam, "slim", sep="_"), rowScale=24, filetype=filetype);
  if (class(resmat) == "matrix") {
    plot_heatmap(resmat, paste(nam, sep="_"), rowScale=24, filetype=filetype);
    plot_heatmap_horiz(resmat, paste(nam, sep="_"), filetype=filetype);
#    if (length(which(colSums(resmat) > 0)) > 0) 
#      plot_heatmap(resmat[,which(colSums(resmat) > 0)], paste(nam, "tall_nonzero", sep="_"), rowScale=24, filetype=filetype);
  }
}

########################################################################################################################

## First letters only
#abbrev <- sapply(names(res_all), function(x) paste(sapply(strsplit(x, " "), substring, 1, 1), collapse="")) 
## Capitals only
#abbrev <- sapply(names(res_all), function(x) paste(sapply(strsplit(x, ""), function(y) y[y %in% LETTERS]), collapse=""));
# Capitals and following lower case letters
abbrev <- sapply(names(res_all), function(x) {
  paste(sapply(strsplit(x, ""), function(y) { 
    sel <- which(y %in% LETTERS);
    sel <- sort(unique(c(sel, sel+1)));
    sel <- sel[y[sel] %in% LETTERS | y[sel] %in% letters]
    y[sel]
  }), collapse="")
})

rownames_filt <- as.character(unlist(sapply(names(res_all), function(x) paste(abbrev[x], rownames(res_all[[x]]), sep=" - "))));

res_all_concat <- do.call(rbind, res_all)
resmat <- as.matrix(res_all_concat, ncol=k);

#resmat_row_ord <- get_optimal_ordering(dist(resmat));
#resmat_col_ord <- get_optimal_ordering(dist(t(resmat)));
resmat_row_ord <- get_optimal_ordering(dist(resmat > 2, method="binary"));
if(exists("resmat_hardcoded", envir=.GlobalEnv)) {
  resmat_col_ord <- resmat_hardcoded 
} else {
  resmat_col_ord <- get_optimal_ordering(dist(t(resmat > 2), method="binary"));
}
resmat <- resmat[resmat_row_ord, resmat_col_ord];

rownames_filt <- rownames_filt[resmat_row_ord];
save(resmat, rownames_filt, file="resmat_rownames_filt.RData")

wordlets <- get_summary_terms(rownames(resmat), 40)
wordlets_pos <- seq(1, nrow(resmat), length.out=length(wordlets)+1);
wordlets_pos <- wordlets_pos[-1] - diff(wordlets_pos)/2;

library(cba)
distmat <- dist(t(resmat > 2), method="binary");
save(distmat, file="distmat_GREAT.RData")
#distmat <- dist(t(resmat))
hc <- hclust(distmat, method="average");
optim <- order.optimal(distmat, hc$merge);
hc$order <- optim$order
hc$merge <- optim$merge;
  
cols=c("white", "yellow", "orange", "red");
cutres <- cutree(hc, k=21)
cutres[cutres %in% which(tabulate(cutres) == 1)] <- 0;
newidx <- as.numeric(factor(as.character(cutres)));
rowcols <- c("white", brewer.pal(max(newidx)-1, "Set1"))[newidx];
  
plotwidth=nrow(resmat)/400;
cexRow <- 1/log10(ncol(resmat)) * 1.5;

for (i in 1:length(abbrev)) {
  print(names(abbrev)[i])
  wordlets_sub <- get_summary_terms(sapply(strsplit(rownames_filt, " - "), function(y) if (y[1] == abbrev[i]) y[2] else NA), 40);

  plotfile(paste(outputdir, "/heatmap_pvals_horiz_all_w_dendro_annot_", abbrev[i], sep=""), height=7, width=plotwidth, type=filetype)
  heatmap.scale(t(resmat), height=8, width=plotwidth, Rowv=as.dendrogram(hc), Colv=NA, 
      col=cols, breaks=c(0,2,3,4,max(c(5, max(resmat[is.finite(resmat)])+1))), main=paste("Annotated using", names(abbrev)[i]),
      scale="none", margins=c(24,8), labCol=NA, cexRow=cexRow, RowSideColors=rowcols, cex.main=0.8,
     add.expr={box(lwd=1); axis(side=1, at=wordlets_pos, labels=wordlets_sub, las=2, cex.axis=1/log10(length(wordlets))*1.5, tick=FALSE, line=-0.5)});
     #add.expr={box(lwd=1); text(x=wordlets_pos, par("usr")[3]-0.5, labels=wordlets, srt=45, pos=2, xpd=TRUE, cex=cexRow)});
  #    add.expr={box(lwd=1); abline(v=cumsum(bnds)+0.5, lwd=1)}));
  #    add.expr={grid(ncol(resmat), nrow(resmat), lty=1, col=1, lwd=0.1); box(lwd=0.1);}));
  dev.off()
} 
  
for (i in 1:length(abbrev)) {
  print(names(abbrev)[i])
  wordlets_sub <- get_summary_terms(sapply(strsplit(rownames_filt, " - "), function(y) if (y[1] == abbrev[i]) y[2] else NA), 40);

  plotfile(paste(outputdir, "/heatmap_pvals_horiz_all_w_annot_", abbrev[i], sep=""), height=7, width=plotwidth, type=filetype)
  heatmap.scale(t(resmat), height=8, width=plotwidth, Rowv=NA, Colv=NA, 
      col=cols, breaks=c(0,2,3,4,max(c(5, max(resmat[is.finite(resmat)])+1))), main=paste("Annotated using", names(abbrev)[i]),
      scale="none", margins=c(24,8), labCol=NA, cexRow=cexRow, cex.main=0.8,
     add.expr={box(lwd=1); axis(side=1, at=wordlets_pos, labels=wordlets_sub, las=2, cex.axis=1/log10(length(wordlets))*1.5, tick=FALSE, line=-0.5)});
     #add.expr={box(lwd=1); text(x=wordlets_pos, par("usr")[3]-0.5, labels=wordlets, srt=45, pos=2, xpd=TRUE, cex=cexRow)});
  #    add.expr={box(lwd=1); abline(v=cumsum(bnds)+0.5, lwd=1)}));
  #    add.expr={grid(ncol(resmat), nrow(resmat), lty=1, col=1, lwd=0.1); box(lwd=0.1);}));
  dev.off()
} 
 


############################################################################################################################################## 
### For chromatin state - based plots WM20160225
load("resmat_rownames_filt_states.RData")
for (i in 1:length(abbrev)) {
  print(names(abbrev)[i])
  wordlets_sub <- get_summary_terms(sapply(strsplit(rownames_filt, " - "), function(y) if (y[1] == abbrev[i]) y[2] else NA), 40);

  plotfile(paste(outputdir, "/heatmap_states_horiz_all_w_annot_", abbrev[i], sep=""), height=7, width=plotwidth, type=filetype)
  heatmap.scale(t(resmat_states), height=8, width=plotwidth, Rowv=NA, Colv=NA, 
      col=cols_norm, breaks=seq(0,numstates)+0.5, main=paste("Annotated using", names(abbrev)[i]),
      scale="none", margins=c(24,8), labCol=NA, cexRow=cexRow, cex.main=0.8,
     add.expr={box(lwd=1); axis(side=1, at=wordlets_pos, labels=wordlets_sub, las=2, cex.axis=1/log10(length(wordlets))*1.5, tick=FALSE, line=-0.5)});
     #add.expr={box(lwd=1); text(x=wordlets_pos, par("usr")[3]-0.5, labels=wordlets, srt=45, pos=2, xpd=TRUE, cex=cexRow)});
  #    add.expr={box(lwd=1); abline(v=cumsum(bnds)+0.5, lwd=1)}));
  #    add.expr={grid(ncol(resmat), nrow(resmat), lty=1, col=1, lwd=0.1); box(lwd=0.1);}));
  dev.off()
} 
  
resmat_states_row_ord <- get_optimal_ordering(dist(resmat_states > 0, method="binary"));
resmat_states_reord <- resmat_states[resmat_states_row_ord, ];
for (i in 1:length(abbrev)) {
  print(names(abbrev)[i])
  wordlets_sub <- get_summary_terms(sapply(strsplit(rownames_filt_states[resmat_states_row_ord], " - "), function(y) if (y[1] == abbrev[i]) y[2] else NA), 40);

  plotfile(paste(outputdir, "/heatmap_states_reord_fixedcols_horiz_all_w_annot_", abbrev[i], sep=""), height=7, width=plotwidth, type=filetype)
  heatmap.scale(t(resmat_states_reord), height=8, width=plotwidth, Rowv=NA, Colv=NA, 
      col=cols_norm, breaks=seq(0,numstates)+0.5, main=paste("Annotated using", names(abbrev)[i]),
      scale="none", margins=c(24,8), labCol=NA, cexRow=cexRow, cex.main=0.8,
     add.expr={box(lwd=1); axis(side=1, at=wordlets_pos, labels=wordlets_sub, las=2, cex.axis=1/log10(length(wordlets))*1.5, tick=FALSE, line=-0.5)});
     #add.expr={box(lwd=1); text(x=wordlets_pos, par("usr")[3]-0.5, labels=wordlets, srt=45, pos=2, xpd=TRUE, cex=cexRow)});
  #    add.expr={box(lwd=1); abline(v=cumsum(bnds)+0.5, lwd=1)}));
  #    add.expr={grid(ncol(resmat), nrow(resmat), lty=1, col=1, lwd=0.1); box(lwd=0.1);}));
  dev.off()
} 
  

whichmax_vals <- apply(apply(resmat_states, 2, tabulate, nbins=numstates), 2, which.max)
max_vals <- apply(apply(resmat_states, 2, tabulate, nbins=numstates), 2, max)
whichmax_vals[which(max_vals == 0)] <- NA
ord <- order(whichmax_vals, -max_vals)
resmat_states_row_ord2 <- order(apply(resmat[,ord], 1, which.max));
resmat_states_reord2 <- resmat_states[resmat_states_row_ord2, rev(ord)];

for (i in 1:length(abbrev)) {
  print(names(abbrev)[i])
  wordlets_sub <- get_summary_terms(sapply(strsplit(rownames_filt_states[resmat_states_row_ord2], " - "), function(y) if (y[1] == abbrev[i]) y[2] else NA), 40);

  plotfile(paste(outputdir, "/heatmap_states_reord2_fixedcols_horiz_all_w_annot_", abbrev[i], sep=""), height=7, width=plotwidth, type=filetype)
  heatmap.scale(t(resmat_states_reord2), height=8, width=plotwidth, Rowv=NA, Colv=NA, 
      col=cols_norm, breaks=seq(0,numstates)+0.5, main=paste("Annotated using", names(abbrev)[i]),
      scale="none", margins=c(24,8), labCol=NA, cexRow=cexRow, cex.main=0.8, RowSideColors=cols_norm[whichmax_vals[rev(ord)]],
     add.expr={box(lwd=1); axis(side=1, at=wordlets_pos, labels=wordlets_sub, las=2, cex.axis=1/log10(length(wordlets))*1.5, tick=FALSE, line=-0.5)});
     #add.expr={box(lwd=1); text(x=wordlets_pos, par("usr")[3]-0.5, labels=wordlets, srt=45, pos=2, xpd=TRUE, cex=cexRow)});
  #    add.expr={box(lwd=1); abline(v=cumsum(bnds)+0.5, lwd=1)}));
  #    add.expr={grid(ncol(resmat), nrow(resmat), lty=1, col=1, lwd=0.1); box(lwd=0.1);}));
  dev.off()
} 
 
for (i in 1:numstates) { 
  sel <- which(resmat_states_reord2 == i, arr.ind=T)
  dat <- resmat_states_reord2[unique(sel[,"row"]), unique(sel[,"col"]),drop=FALSE]

  if (ncol(dat) < 2 | nrow(dat) < 2) next;

  wordlets <- get_summary_terms(rownames(dat), ceiling(nrow(dat)/100), filtering=FALSE)
  wordlets_pos <- seq(1, nrow(dat), length.out=length(wordlets)+1);
  wordlets_pos <- wordlets_pos[-1] - diff(wordlets_pos)/2;

  plotfile(paste(outputdir, "/heatmap_states_reord2_fixedcols_horiz_all_state_", i, sep=""), height=7, width=10, type=filetype)
  heatmap.scale(t(dat), height=8, width=10, Rowv=NA, Colv=NA, 
      col=cols_norm, breaks=seq(0,numstates)+0.5, main=paste("State", state_names[i]),
      scale="none", margins=c(24,8), labCol=NA, cexRow={1/log10(ncol(dat)) * 1.5}, cex.main=0.8,
     add.expr={box(lwd=1); axis(side=1, at=wordlets_pos, labels=wordlets, las=2, cex.axis=1/log10(length(wordlets))*1.5, tick=FALSE, line=-0.5)});
     #add.expr={box(lwd=1); text(x=wordlets_pos, par("usr")[3]-0.5, labels=wordlets, srt=45, pos=2, xpd=TRUE, cex=cexRow)});
  #    add.expr={box(lwd=1); abline(v=cumsum(bnds)+0.5, lwd=1)}));
  #    add.expr={grid(ncol(resmat), nrow(resmat), lty=1, col=1, lwd=0.1); box(lwd=0.1);}));
  dev.off()
} 

load("Func_Factors.RData")
dat_factors <- as.matrix(d[,4:ncol(d)])
rownames(dat_factors) <- d$func;
for (i in 1:numstates) { 
  sel <- which(resmat_states_reord2 == i, arr.ind=T)
  dat <- resmat_states_reord2[unique(sel[,"row"]), unique(sel[,"col"]),drop=FALSE]

  if (ncol(dat) < 2 | nrow(dat) < 2) next;

  dat_factors_state <- round(dat_factors[colnames(dat),], 1)
  dat_factors_state <- dat_factors_state[,colSums(dat_factors_state) != 0]

  whichmax_vals <- apply(dat_factors_state, 2, which.max)
  max_vals <- apply(dat_factors_state, 2, max)
  ord <- order(whichmax_vals, -max_vals)

  plotfile(paste(outputdir, "/heatmap_states_reord2_fixedcols_horiz_CRs_state_", i, sep=""), height=7, width=3, type=filetype)
  heatmap.scale(dat_factors_state[,ord], col=colorpanel(100, low="white", high="black"), Rowv=NA, Colv=NA, scale="none", margins=c(10,8), height=7, width=3);
  dev.off()
} 


resmat_states_reord3 <- resmat_states_reord2[rowSums(resmat_states_reord2>0)>1, colSums(resmat_states_reord2>0)>1]
neword <- order(sapply(apply(resmat_states_reord3, 1, tabulate), which.max))
neword2 <- rev(order(sapply(apply(resmat_states_reord3, 2, tabulate), which.max)))
wordlets <- get_summary_terms(rownames(resmat_states_reord3[neword,neword2]), ceiling(nrow(resmat_states_reord3[neword,neword2])/100), filtering=FALSE)
wordlets_pos <- seq(1, nrow(resmat_states_reord3[neword,neword2]), length.out=length(wordlets)+1);
wordlets_pos <- wordlets_pos[-1] - diff(wordlets_pos)/2;
plotfile(paste(outputdir, "/heatmap_states_reord3", sep=""), height=7, width=13, type=filetype)
heatmap.scale(t(resmat_states_reord3[neword,neword2]), col=cols_norm, breaks=seq(0,numstates)+0.5, Rowv=NA, Colv=NA, height=9, width=15,
        scale="none", margins=c(24,8), labCol=NA, cexRow={1/log10(ncol(resmat_states_reord3[neword,neword2])) * 1.5}, cex.main=0.8,
        RowSideColors=rev(cols_norm[sort(sapply(apply(resmat_states_reord3, 2, tabulate), which.max))]),
        ColSideColors=cols_norm[sort(sapply(apply(resmat_states_reord3, 1, tabulate), which.max))],
        add.expr={box(lwd=1); axis(side=1, at=wordlets_pos, labels=wordlets, las=2, cex.axis=1/log10(length(wordlets))*1.5, tick=FALSE, line=-0.5)});
dev.off()

dat_factors3 <- dat_factors[colnames(resmat_states_reord3)[neword2],]
#dat_factors3 <- dat_factors[colSums(resmat_states_reord2>0)>1,];
dat_factors3 <- dat_factors3[,colSums(dat_factors3 > 0.5) > 0];

plotfile(paste(outputdir, "/heatmap_factors_reord3", sep=""), height=8, width=7, type=filetype)
heatmap.scale(dat_factors3[,rev(order(apply(dat_factors3, 2, which.max), apply(dat_factors3, 2, max)))], 
        col=colorpanel(100, low="white", high="black"), Rowv=NA, Colv=NA, height=9, width=7,
        scale="none", margins=c(24,8), cex.main=0.8,
        add.expr={box(lwd=1); grid(nx=ncol(dat_factors3), ny=nrow(dat_factors3), lwd=0.5,lty=1,col="grey")})
dev.off()

######## Start over.
### Filter resmat to be highly significant (-log10(p) > 6).
sel_strict_rows <- which(rowSums(resmat>6)>1);
sel_strict_cols <- which(colSums(resmat[sel_strict_rows,]>6)>1);
#sel_strict_rows <- sel_strict_rows[which(rowSums(resmat_states[sel_strict_rows, sel_strict_cols]>0)>1)]
#sel_strict_cols <- sel_strict_cols[which(colSums(resmat_states[sel_strict_rows, sel_strict_cols]>0)>1)]
resmat_states_strict <- resmat_states[sel_strict_rows,sel_strict_cols]
rownames_filt_states_strict <- rownames_filt_states[sel_strict_rows]

ord_strict_rows <- order(sapply(apply(resmat_states_strict, 1, tabulate), which.max))
ord_strict_cols <- rev(order(sapply(apply(resmat_states_strict, 2, tabulate), which.max)))
resmat_states_strict_ord <- resmat_states_strict[ord_strict_rows,ord_strict_cols]
wordlets <- get_summary_terms(rownames(resmat_states_strict_ord), ceiling(nrow(resmat_states_strict_ord)/100), filtering=FALSE)
wordlets_pos <- seq(1, nrow(resmat_states_strict_ord), length.out=length(wordlets)+1);
wordlets_pos <- wordlets_pos[-1] - diff(wordlets_pos)/2;
plotfile(paste(outputdir, "/heatmap_states_strict_ord", sep=""), height=7, width=13, type=filetype)
heatmap.scale(t(resmat_states_strict_ord), col=cols_norm, breaks=seq(0,numstates)+0.5, Rowv=NA, Colv=NA, height=9, width=15,
        scale="none", margins=c(24,8), labCol=NA, cexRow={1/log10(ncol(resmat_states_strict_ord)) * 1.5}, cex.main=0.8,
        RowSideColors=cols_norm[sapply(apply(resmat_states_strict_ord, 2, tabulate), which.max)],
        ColSideColors=cols_norm[sapply(apply(resmat_states_strict_ord, 1, tabulate), which.max)],
        add.expr={box(lwd=1); axis(side=1, at=wordlets_pos, labels=wordlets, las=2, cex.axis=1/log10(length(wordlets))*1.5, tick=FALSE, line=-0.5)});
dev.off()

dat_factors_strict <- dat_factors[colnames(resmat_states_strict_ord),]
#dat_factors3 <- dat_factors[colSums(resmat_states_reord2>0)>1,];
dat_factors_strict <- dat_factors_strict[,colSums(dat_factors_strict > 0.5) > 0];

plotfile(paste(outputdir, "/heatmap_factors_strict", sep=""), height=8, width=7, type=filetype)
heatmap.scale(dat_factors_strict[,rev(order(apply(dat_factors_strict, 2, which.max), apply(dat_factors_strict, 2, max)))],
        col=colorpanel(100, low="white", high="black"), Rowv=NA, Colv=NA, height=9, width=7,
        scale="none", margins=c(24,8), cex.main=0.8,
        add.expr={box(lwd=1); grid(nx=ncol(dat_factors_strict), ny=nrow(dat_factors_strict), lwd=0.5,lty=1,col="grey")})
dev.off()




### Filter resmat to be highly significant (-log10(p) > 6).
### ONLY GO Biological Processes!
sel_strict_rows <- which(sapply(strsplit(rownames_filt, " - "), "[[", 1) == "GOBiPr");
sel_strict_rows <- sel_strict_rows[which(rowSums(resmat[sel_strict_rows,]>6)>1)];
sel_strict_cols <- which(colSums(resmat[sel_strict_rows,]>6)>1);
sel_strict_rows <- sel_strict_rows[which(rowSums(resmat_states[sel_strict_rows, sel_strict_cols]>0)>0)]
sel_strict_cols <- sel_strict_cols[which(colSums(resmat_states[sel_strict_rows, sel_strict_cols]>0)>0)]
#sel_strict_rows <- sel_strict_rows[which(rowSums(resmat_states[sel_strict_rows, sel_strict_cols]>0)>1)]
#sel_strict_cols <- sel_strict_cols[which(colSums(resmat_states[sel_strict_rows, sel_strict_cols]>0)>1)]
resmat_states_strict <- resmat_states[sel_strict_rows,sel_strict_cols]
rownames_filt_states_strict <- rownames_filt_states[sel_strict_rows]

ord_strict_rows <- order(sapply(apply(resmat_states_strict, 1, tabulate), which.max))
ord_strict_cols <- rev(order(sapply(apply(resmat_states_strict, 2, tabulate), which.max)))
resmat_states_strict_ord <- resmat_states_strict[ord_strict_rows,ord_strict_cols]
#wordlets <- get_summary_terms(rownames(resmat_states_strict_ord), ceiling(nrow(resmat_states_strict_ord)/100), filtering=FALSE)
wordlets <- get_summary_terms(rownames(resmat_states_strict_ord), 2, filtering=FALSE)
wordlets_pos <- seq(1, nrow(resmat_states_strict_ord), length.out=length(wordlets)+1);
wordlets_pos <- wordlets_pos[-1] - diff(wordlets_pos)/2;
Rowcols <- sapply(apply(resmat_states_strict_ord, 2, tabulate), which.max)
Colcols <- sapply(apply(resmat_states_strict_ord, 1, tabulate), which.max)

plotfile(paste(outputdir, "/heatmap_states_strict_ord_GOBiPr", sep=""), height=7, width=16, type=filetype)
heatmap.scale(t(resmat_states_strict_ord), col=cols_norm, breaks=seq(0,numstates)+0.5, Rowv=NA, Colv=NA, height=9, width=22,
        scale="none", margins=c(24,4), labCol=NA, cexRow={1/log10(ncol(resmat_states_strict_ord)) * 1.5}, cex.main=0.8,
        RowSideColors=cols_norm[Rowcols], ColSideColors=cols_norm[Colcols],
        add.expr={box(lwd=1); axis(side=1, at=wordlets_pos, labels=wordlets, las=2, cex.axis=1/log10(length(wordlets))*1.5, tick=FALSE, line=-0.5); abline(h=which(diff(Rowcols) != 0)+0.5); abline(v=which(diff(Colcols) != 0)+0.5)});
dev.off()

dat_factors_strict <- dat_factors[colnames(resmat_states_strict_ord),]
#dat_factors3 <- dat_factors[colSums(resmat_states_reord2>0)>1,];
dat_factors_strict <- dat_factors_strict[,colSums(dat_factors_strict > 0.5) > 0];

plotfile(paste(outputdir, "/heatmap_factors_strict_GOBiPr", sep=""), height=8, width=7, type=filetype)
heatmap.scale(dat_factors_strict[,rev(order(apply(dat_factors_strict, 2, which.max), apply(dat_factors_strict, 2, max)))],
        col=colorpanel(100, low="white", high="black"), Rowv=NA, Colv=NA, height=9, width=7,
        scale="none", margins=c(24,8), cex.main=0.8,
        add.expr={box(lwd=1); grid(nx=ncol(dat_factors_strict), ny=nrow(dat_factors_strict), lwd=0.5,lty=1,col="grey")})
dev.off()





