#!/usr/bin/R
# Plot available and imputed data for chromImpute:
# 1. Plot sample mark and % missing for all 
# 2. Plot imputation choices (cutoffs) for table to choose which epitopes to impute.
# 3. Plot similarity on top of available data.
# 4. Plot success of imputation of available data (VALIDATION).
# 5. Plot success of imputation jointly with similarity.
# 6. TBD Other plots from original CI paper.
source('~/EPIMAP_ANALYSIS/bin/general_EPIMAP_ANALYSIS.R')
today <- Sys.Date()
today <- format(today, format="%m%d%Y")

# All processing files:
tab <- read.delim('ChromImpute/sample_mark_table.tsv', header=F)
mapping <- read.delim('Annotation/all_submitted_released_biosample_mapping.tsv')
names(tab) <- c('CellType', 'Epitope', 'File')
names(mapping) <- c('CellTypeName','CellType')
tab$Number <- 1
tab <- merge(tab, mapping)

# Reordering by availability: 
cells <- aggregate(File ~ CellType, tab, length)
epitopes <- aggregate(File ~ Epitope, tab, length)
ordCT <- as.character(with(cells,CellType[order(File, decreasing=TRUE)]))
ordE <- as.character(with(epitopes,Epitope[order(File, decreasing=TRUE)]))

# Long to wide table:
tmp <- dcast(tab, Epitope ~ CellType, value.var='Number')
rownames(tmp) <- tmp$Epitope 
wide <- tmp[ordE,ordCT]

# ==============================================
# 1. Plot of available sample/mark combinations:
# ==============================================
wm <- as.matrix(wide)
matcol=c('white','darkblue')
wm[is.na(wm)] <- 0
cols <- (seq(0,nrow(wm),5) - .5)/ (nrow(wm) -1)
rows <- (seq(0,ncol(wm),50) - .5)/ (ncol(wm) -1)

png(paste0(img,'/sample_mark_table_', today, '.png'),res=450,units='in',width=8.5,height=11)
par(mar=c(7,4,3,2))
image(wm, axes=FALSE,main='',col=matcol)
mtext('ENCODE Sample/Mark Combinations', side=3, cex=2.25, line=.75)
mtext(paste0('Biosamples (', ncol(wm), ')'), side=2, cex=2, line=1)
mtext(paste0('ChIP-seq Epitopes and DNase-seq (', nrow(wm), ')'), side=1, cex=2, line=5.5)
grid(nx=nrow(wm), ny=ncol(wm)/5,col='grey',lty='solid',lwd=.25)
text(x=seq(0,1,length.out=nrow(wm)), y=par()$usr[3]-0.01*(par()$usr[4]-par()$usr[3]), labels=rownames(wm), srt=90, adj=1, xpd=TRUE,cex=.75)
abline(h=rows,col='black',lty='solid',lwd=.5)
abline(v=cols,col='black',lty='dashed',lwd=.5)
dev.off()

# ==================================================================================
# 2. Plot imputation choices (cutoffs) for table to choose which epitopes to impute.
# ==================================================================================
N.EPI <- length(ordE)
pct <- rep(0, N.EPI)
cuts <- 1:N.EPI
for (cutoff in cuts){
    mat <- wm[ordE[1:cutoff],]
    if (cutoff == 1){
        pct[cutoff] = sum(mat) / length(mat)
    } else {
        pct[cutoff] = sum(mat) / prod(dim(mat))
    }
}

breaks = c(6,8,13,25,36)
bcols <- (breaks - .5)/ (nrow(wm) -1)
bpct = pct[breaks]
ll <- diff(c(0,breaks))
mids <- ll / 2 + c(0, breaks)[-(length(breaks) +1)]
textpos <- (mids -.5) / (nrow(wm) - 1)
ypos = .55
# ypos = rep(seq(.85,.5,length.out=breaks), length(breaks))[1:length(breaks)]
projtext = paste0(breaks, ' marks (', round(bpct*100, 1), '%)')

png(paste0(img,'/sample_mark_table_cutoffs_', today,'.png'),res=450,units='in',width=8.5,height=11)
par(mar=c(7,4,3,2))
image(wm, axes=FALSE,main='',col=matcol)
mtext('ENCODE Sample/Mark Combinations', side=3, cex=2.25, line=.75)
mtext(paste0('Biosamples (', ncol(wm), ')'), side=2, cex=2, line=1)
mtext(paste0('ChIP-seq Epitopes and DNase-seq (', nrow(wm), ')'), side=1, cex=2, line=5.5)
grid(nx=nrow(wm), ny=ncol(wm)/5,col='grey',lty='solid',lwd=.25)
text(x=seq(0,1,length.out=nrow(wm)), y=par()$usr[3]-0.01*(par()$usr[4]-par()$usr[3]), labels=rownames(wm), srt=90, adj=1, xpd=TRUE,cex=.75)
abline(h=rows,col='black',lty='solid',lwd=.5)
abline(v=bcols,col='darkred',lty='dashed',lwd=2)
text(x=textpos, y=ypos, labels=projtext, srt=90, adj=0, xpd=TRUE, col=alpha('darkred',.85), cex=2.5)
dev.off()

# ============================================
# 3. Plot similarity on top of available data.
# 4. Plot success of imputation of available data (VALIDATION).
# 5. Plot success of imputation jointly with similarity.
# ============================================

