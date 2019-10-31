#!/usr/bin/R
# Return distances from a 
source('~/EPIMAP_ANALYSIS/bin/general_EPIMAP_ANALYSIS.R')
args=(commandArgs(TRUE))
file = args[1]
project = args[2]

# Read in data: 
cont <- read.table(file,header=F)
names(cont) <- c('pref','comp','x00','x01','x10','x11', 'NR')
prefixes = sort(unique(c(as.character(cont$pref),
                         as.character(cont$comp))))
cont$pref <- factor(cont$pref, levels=prefixes)
cont$comp <- factor(cont$comp, levels=prefixes)

celltypes <- sort(unique(matrix(unlist(strsplit(prefixes,"/")),nrow=2)[2,]))

# Marginals:
cont$total <- cont$x00 + cont$x01 + cont$x10 + cont$x11
cont$p0 <- (cont$x00 + cont$x01)/cont$total
cont$p1 <- (cont$x10 + cont$x11)/cont$total
cont$c0 <- (cont$x00 + cont$x10)/cont$total
cont$c1 <- (cont$x01 + cont$x11)/cont$total

# --------------------------------------------
# Similarity and information content measures:
# --------------------------------------------
# 1. Jaccard: 
cont$jaccard <- with(cont, x11 / (x11 + x01 + x10))

# 2. Mutual information:
# -- sum p(x,y) log(p(x,y)/(p(x)p(y))) for (x,y) in (0,0), (0,1), (1,0), (1,1).
cont$mutual_inf <- with(cont,
                        x00/total * log((x00/total) / (p0 * c0)) + 
                            x10/total * log((x10/total) / (p1 * c0)) + 
                            x01/total * log((x01/total) / (p0 * c1)) + 
                            x11/total * log((x11/total) / (p1 * c1)))

# NOTE: Entropy measures are not particularly useful.
# Entropy H(X) - for each prefix and comp
cont$p.entropy <- with(cont, - log(p0) * p0 - log(p1) * p1)
cont$c.entropy <- with(cont, - log(c0) * c0 - log(c1) * c1)

# Conditional entropy: (non-symmetric)
cont$cond.entropy <- with(cont,
                          x00/total * log(p0/(x00/total)) + 
                              x10/total * log(p1/(x10/total)) + 
                              x01/total * log(p0/(x01/total)) + 
                              x11/total * log(p1/(x11/total)))

# Conditional entropy: (non-symmetric)
cont$cond.entropy2 <- with(cont,
                           x00/total * log(c0/(x00/total)) + 
                               x10/total * log(c0/(x10/total)) + 
                               x01/total * log(c1/(x01/total)) + 
                               x11/total * log(c1/(x11/total)))


# ---------------------------------
# Matrix format for these measures:
# With prefix X comp
# ---------------------------------
N = length(prefixes)
dfid = cbind(cont$pref,cont$comp)
rev.dfid = cbind(cont$comp,cont$pref) # only for symmetric measures
# 1. Jaccard:
jaccard <- matrix(0, nrow=N, ncol=N, dimnames=list(prefixes, prefixes))
jaccard[dfid] <- cont$jaccard
jaccard[rev.dfid] <- cont$jaccard
# Reorder:
jaccard <- reord(jaccard) 
jaccard <- jaccard[,rownames(jaccard)]

# Image: 
png(paste0(img,'/jaccard_distance_cont_',project,'.png'),res=450,units='in',width=12.5,height=11)
par(mar=c(2,14,2,2))
image(jaccard, axes=FALSE, main="All Datasets (Jaccard Similarity)",col=colramp)
grid(nx=N, ny=N, col='grey',lty='solid',lwd=.25)
text(y=seq(0,1,length.out=N), x=par()$usr[3]-0.001*(par()$usr[4]-par()$usr[3]), labels=colnames(jaccard), srt=0, adj=1, xpd=TRUE,cex=0.5)
dev.off()

write.table(jaccard, paste0('jaccard_distance_cont_',project,'.tsv'), sep="\t", quote=F)

# 2. Mutual Inf:
# --------------
mutual_inf <- matrix(0, nrow=N, ncol=N, dimnames=list(prefixes, prefixes))
mutual_inf[dfid] <- cont$mutual_inf
mutual_inf[rev.dfid] <- cont$mutual_inf
# Reorder:
mutual_inf <- reord(mutual_inf) 
mutual_inf <- mutual_inf[,rownames(mutual_inf)]

# Image: 
png(paste0(img,'/mutual_inf_distance_cont_',project,'.png'),res=450,units='in',width=12.5,height=11)
par(mar=c(2,14,2,2))
image(mutual_inf, axes=FALSE, main="All Datasets (Mutual Information)",col=colramp)
grid(nx=N, ny=N, col='grey',lty='solid',lwd=.25)
text(y=seq(0,1,length.out=N), x=par()$usr[3]-0.001*(par()$usr[4]-par()$usr[3]), labels=colnames(mutual_inf), srt=0, adj=1, xpd=TRUE,cex=0.5)
dev.off()

write.table(mutual_inf, paste0('mutual_inf_distance_cont_',project,'.tsv'), sep="\t", quote=F)

# -------------
# By cell type:
# -------------
for (ct in celltypes){
    print(ct)
    ids <- grep(paste0("/",ct), colnames(jaccard))

    png(paste0(img,'/jaccard_distance_cont_',project,'_',ct,'.png'),res=450,units='in',width=12.5,height=11)
    par(mar=c(2,14,2,2))
    image(jaccard[ids,ids], axes=FALSE, main=paste(ct,"(Jaccard Similarity)"),col=colramp)
    grid(nx=length(ids), ny=length(ids), col='grey',lty='solid',lwd=.25)
    text(y=seq(0,1,length.out=length(ids)), x=par()$usr[3]-0.01*(par()$usr[4]-par()$usr[3]), labels=colnames(jaccard)[ids], srt=0, adj=1, xpd=TRUE,cex=0.75)
    dev.off()

    png(paste0(img,'/mutual_inf_distance_cont_',project,'_',ct,'.png'),res=450,units='in',width=12.5,height=11)
    par(mar=c(2,14,2,2))
    image(mutual_inf[ids,ids], axes=FALSE, main=paste(ct,"(Mutual Information)"),col=colramp)
    grid(nx=length(ids), ny=length(ids), col='grey',lty='solid',lwd=.25)
    text(y=seq(0,1,length.out=length(ids)), x=par()$usr[3]-0.01*(par()$usr[4]-par()$usr[3]), labels=colnames(mutual_inf)[ids], srt=0, adj=1, xpd=TRUE,cex=0.75)
    dev.off()

}

write.table(cont, paste0('full_distance_cont_save_',project,'.tsv'), sep="\t", quote=F)

