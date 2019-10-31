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
library(colormap)
source(paste0(bindir,'general_EPIMAP_ANALYSIS.R'))
source(paste0(bindir, 'great_enrichment/general_WM_functions.R'))
source(paste0(bindir, 'great_enrichment/auxfunctions_GREAT.R'))

# --------------------------------------------
# NOTE: Starts with original plotting scripts.
# Requires arguments for the following code:
# --------------------------------------------
source(paste0(bindir, 'plot_GREAT_results.R'))


# ----------------------------------------------
### For chromatin state - based plots WM20160225
# ----------------------------------------------
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





