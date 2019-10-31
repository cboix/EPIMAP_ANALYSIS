source("../../code/general.R")
library(ade4)
library(gplots)
library(RColorBrewer)

get_summary_terms <- function(dat, binsize=50) {
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
    ranks[tab < 2] <- NA;

    term_score <- sapply(term_words[i:j], function(x) {
      mean(ranks[names(ranks) %in% x], na.rm=T);
    });

    val <- which.max(term_score);
    terms_sel <- c(terms_sel, ifelse(length(val)==0, "", dat[i:j][val]));
  }

  invisible(terms_sel);
}

load("resmat_rownames_filt.RData")# resmat & rownames_filt
dat_genesets <- t(resmat);
#GOBiPr_sel <- which(sapply(sapply(rownames_filt, strsplit, " - "), "[[", 1) == "GOBiPr"); # GO Biological Processes only!
#dat_genesets <- dat_genesets[,GOBiPr_sel]
#rownames_filt <- rownames_filt[GOBiPr_sel];

cnames <- colnames(dat_genesets)
dup_hits <- which(duplicated(cnames))
cnames[dup_hits] <- paste(cnames[dup_hits], "2");
colnames(dat_genesets) <- cnames;
cnames <- colnames(dat_genesets)
dup_hits <- which(duplicated(cnames))
cnames[dup_hits] <- paste(cnames[dup_hits], "2");
colnames(dat_genesets) <- cnames;

load("Func_Factors.RData")
dat_factors <- as.matrix(d[,4:ncol(d)])
rownames(dat_factors) <- d$func;
factor_ord <- c("HDAC6.1", "ZNF263", "BRF2", "ILF2", "PHF8.1", "GABPB1", "YY1", "IRF1", "CCNE1", "CCNT2", "CHD2", "E2F4", "E2F6", "HMGN3", "MAX", "MAZ", 
  "MXI1", "MYC", "NRF1", "SP1", "SP2", "THAP1", "XRCC4", "BCL3", "BCLAF1", "BRF1", "CHD1.1", "ELK1", "ETS1", "GTF2F1", "H3K4me3", "HDAC1.1", "HDAC8", 
  "KDM5B", "NR2C2", "PML", "POLR2A", "RbBP5", "RBBP5", "RDBP", "SAP30.1", "SIN3A", "STAT1", "STAT2", "TAF1", "TAF7", "TBL1XR1", "UBTF", "ARID3A", 
  "EP300", "FOSL1", "GATA1", "JUN", "SIRT6.1", "SMARCB1", "STAT5A", "ATF3", "FOS", "NFE2", "NFYA", "SIX5", "USF1", "USF2", "NFYB", "HDAC2.1", 
  "SMARCA4", "RFX5", "KDM5A_JARID1A_RbBP2", "WDR77", "MLL2", "NCOR1", "KDM3A_JMJD1A_JHDM2A", "CBX7", "PHF8", "SMARCA4_BRG1", "SAP30", "CHD7", 
  "CHD1", "KAT2B_PCAF", "SIRT6", "WDR5", "H3K9me3", "CBX3_HP1.gamma", "ZNF274", "TRIM28", "CBX3", "SETDB1", "ZBTB33", "PAF1", "CHD4_Mi2", "CARM1_PRMT4", 
  "ESR1", "RNApol2S5p1", "NR4A1", "H3K36me3", "WHSC1_NSD2", "H3K27ac", "H3K79me2", "H3K4me2", "H3K4me1", "H2A.Z", "LSD1_BHC110", "KDM6A_UTX", 
  "SMARCA5_SNF2H_ISWI", "KDM5B_JARID1B_PLU1", "ASH2L_ASH2", "BACH1", "MAFF", "MAFK", "EGR1", "NR2F2", "GATA2", "TAL1", "TEAD4", "DNase", "ATF1", 
  "JUNB", "JUND", "CEBPB", "CTCF", "CTCF.1", "CTCFL", "RAD21", "SMC3", "ZNF143", "BHLHE40", "SPI1", "ELF1", "ZBTB7A", "SRF", "BDP1", "GTF2B", 
  "MEF2A", "POLR3A", "POLR3G", "TBP", "BRD4", "RING1", "SIRT2", "SUZ12", "H3K27me3", "EED", "EZH2", "EZH2.1", "MLL5", "CBX8", "CBX2", "RNF2_RING1B", 
  "H3K9me2", "EP300_KAT3B_p300", "JMJD1C", "CREBBP_KAT3A_CBP", "SIRT1", "CtBP1", "HDAC6", "HDAC2", "HDAC1", "SETDB1_ESET", "EHMT2_G9A", "RCOR1", 
  "REST.1", "REST", "KDM5C_JARID1C_SMCX") 
dat_factors <- dat_factors[rownames(dat_genesets), factor_ord];

pca_factors <- dudi.pca(dat_factors, scale = TRUE, scan = FALSE)
pca_genesets <- dudi.pca(dat_genesets, scale = TRUE, scan = FALSE)

coin1 <- coinertia(pca_factors, pca_genesets, scan = FALSE, nf = 2)
names(coin1)

tab <- coin1$tab;
tab_row_ord <- get_optimal_ordering(dist(tab));
tab_col_ord <- get_optimal_ordering(dist(t(tab)));

plotfile("test", width=32, type="pdf")
heatmap.scale(t(as.matrix(tab[tab_row_ord,tab_col_ord])), scale="none", cexCol=0.2, cexRow=0.2, 
  width=32, col=colorpanel(100, low="white", high="black"))
dev.off()

plotwidth=max(10, nrow(tab)/400);
cexRow <- 1/log10(ncol(tab));
filetype="pdf"

hc <- hclust(dist(t(tab)), method="average");
optim <- order.optimal(dist(t(tab)), hc$merge);
hc$order <- optim$order
hc$merge <- optim$merge;

cutres <- cutree(hc, k=18)
cutres[cutres %in% which(tabulate(cutres) == 1)] <- 0;
newidx <- as.numeric(factor(as.character(cutres)));
rowcols <- c("white", brewer.pal(max(newidx)-1, "Set1"))[newidx];

wordlets <- get_summary_terms(rownames(tab)[tab_row_ord], min(plotwidth,20))
wordlets_pos <- seq(1, nrow(tab), length.out=length(wordlets)+1);
wordlets_pos <- wordlets_pos[-1] - diff(wordlets_pos)/2;
wordlets_sub <- get_summary_terms(sapply(strsplit(rownames_filt[tab_row_ord], " - "), function(y) if (y[1] == "GOBiPr") y[2] else NA), min(plotwidth,20));
#wordlets_sub <- get_summary_terms(sapply(strsplit(rownames_filt[tab_row_ord], " - "), function(y) if (y[1] == "GOBiPr" | y[1] == "GOMoFu" | y[1] == "GOCeCo") y[2] else NA), 20);

outputdir <- ".";
plotfile(paste(outputdir, "/heatmap_coin1", sep=""), height=7, width=plotwidth, type=filetype)
heatmap.scale(t(tab[tab_row_ord,]), height=8, width=plotwidth, Colv=NA, Rowv=as.dendrogram(hc), RowSideColors=rowcols,
   col=colorpanel(100, low="white", high="black"), scale="none", margins=c(20,8), labCol=NA, cexRow=cexRow,
   add.expr={box(lwd=1); axis(side=1, at=wordlets_pos, labels=wordlets_sub, las=2, cex.axis=1/log10(length(wordlets))*1.5, tick=FALSE, line=-0.5)});
dev.off()

outputdir <- ".";
plotfile(paste(outputdir, "/heatmap_coin1_all_wordlets", sep=""), height=7, width=plotwidth, type=filetype)
heatmap.scale(t(tab[tab_row_ord,]), height=8, width=plotwidth, Colv=NA, Rowv=as.dendrogram(hc), RowSideColors=rowcols,
   col=colorpanel(100, low="white", high="black"), scale="none", margins=c(20,8), labCol=NA, cexRow=cexRow,
   add.expr={box(lwd=1); axis(side=1, at=wordlets_pos, labels=wordlets, las=2, cex.axis=1/log10(length(wordlets))*1.5, tick=FALSE, line=-0.5)});
dev.off()





