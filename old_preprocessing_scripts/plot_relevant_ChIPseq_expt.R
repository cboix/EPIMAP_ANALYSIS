#!/usr/bin/R
# Plot the processing status report for ENCODE3/ROADMAP DATA
source('~/EPIMAP_ANALYSIS/bin/general_EPIMAP_ANALYSIS.R')

read.table('relevant_experiments',header=F) -> df
colnames(df) <- c('EPITOPE','CT')
as.character(sort(unique(df$EPITOPE))) -> epitopes
as.character(sort(unique(df$CT))) -> ct

NE=length(epitopes)
NC=length(ct)
mat <- matrix(0, nrow=NE, ncol=NC, dimnames=list(epitopes,ct))

df$EPITOPE <- factor(df$EPITOPE, levels=epitopes)
df$CT <- factor(df$CT, levels=ct)
dfid <- as.matrix(df)
mat[dfid] <- 1

cols = 4
cuts = round(seq(0, NE, length.out=cols+1))
rangelist = list()
matlist = list()
for (i in 1:cols){ 
    rangelist[[i]] = (cuts[i]+1):cuts[i+1]
    matlist[[i]] = mat[rangelist[[i]],]
}


png(paste0(img,'/relevant_ChIP-seq_experiments.png'),res=450,units='in',height=17,width=11)
layout(matrix(1:cols, 1, cols), widths=rep(1, cols), TRUE)
par(yaxs="i")
par(xaxs="i")
for (i in 1:cols){ 
    wm = t(matlist[[i]])
    par(mar=c(5.5,8,1,0.5))
    image(wm, axes=FALSE, main="", col=c('white','royalblue'))
    grid(nx=nrow(wm), ny=ncol(wm), col='grey', lty='solid', lwd=.35)
    text(x=seq(0,1,length.out=nrow(wm)), y=par()$usr[3]-0.01*(par()$usr[4]-par()$usr[3]), labels=rownames(wm), srt=90, adj=1, xpd=TRUE,cex=1)
text(y=seq(0,1,length.out=ncol(wm)), x=par()$usr[3]-0.15*(par()$usr[4]-par()$usr[3]), labels=colnames(wm), srt=0, adj=1, xpd=TRUE,cex=.75)
}
dev.off()

