#!/usr/bin/R
# -----------------------------------------------
# Given intersection from BW and PK file:
# 1. Collect aggregate statistics
# 2. Plot each epigenome
# 3. Save as matrix
# -----------------------------------------------

domain = system("hostname -d", intern=TRUE)
if (domain == 'broadinstitute.org'){
    bindir='~/data/EPIMAP_ANALYSIS/bin/'
} else {
    bindir='~/EPIMAP_ANALYSIS/bin/'
}
source(paste0(bindir, 'general_EPIMAP_ANALYSIS.R'))
source(paste0(bindir, 'auxiliary_chromImpute_functions.R'))
source(paste0(bindir, 'load_metadata.R'))

options(scipen=45) # so we dont get issues writing integers into bedfiles
today <- format(Sys.time(), "%m%d%y")

# Args
# filename='BSS00001_pm250.tsv.gz'
# tagline='BSS00001_pm250'
# statsfile='BSS00001_pm250_stats.tsv'
# plotfile='BSS00001_pm250_heatmap.png'

# Arguments:
args=(commandArgs(TRUE))
if (length(args)==0) {
    stop("No arguments supplied. Need:
          filename tagline statsfile plotfile")
} else {        
    filename = args[1]
    tagline = args[2]
    statsfile = args[3]
    plotfile = args[4]
}

# Read in intersection:
df = read.delim(gzfile(filename), header=F, sep="\t")
names(df) = c('id','start','end','value')

# Make matrix:
# NOTE: only works because of these midpts
df$bin = round((df$end + 238)/25 + 1)
wide = spread(df[,c('id','bin','value')], bin, value)
mat = as.matrix(wide[,-1])
rownames(mat) = wide[,1]

# Get volume of center:
marg = apply(mat[,8:14], 1, sum) + 0.5 * apply(mat[,c(7,15)], 1, sum)
marg = marg / 200 * 25
mdf = data.frame(id=rownames(mat), marg=marg)

# Order matrix by volume of center
ord = order(marg, decreasing=FALSE)
# remat = log2(mat[ord,])
remat = mat[ord,]
remat[remat > 8] = 8

frac2 = sum(marg > 2) / length(marg)
label = sprintf("%0.1f%%\n\n%0.0d\nout of\n%0.0d",
                frac2 * 100, sum(marg > 2), length(marg))

# ------------------------------
# Output stats (before we plot):
# ------------------------------
stats = c('tot.dhs' = length(marg), 
          'frac2' = round(100*frac2, 2),
          'gt2' = sum(marg > 2),
          'gt3' = sum(marg > 3),
          'gt4' = sum(marg > 4),
          'gt5' = sum(marg > 5))

statsdf = data.frame(metric = names(stats), 
                     value = stats, 
                     id = tagline)

write.table(statsdf, statsfile, sep="\t", quote=F, row.names=F)

# ------------------------------------
# Save ordered matrix: (for plotting?)
# ------------------------------------
# matfile = sub(".tsv", "_mat.tsv", filename)
# write.table(remat, gzfile(matfile), sep="\t", quote=F)

# ------------
# Plot matrix:
# ------------
png(plotfile, width=1, height=7, res=150, units='in')
par(mar=c(1.5,.25, 1.5, .25))
par(yaxs="i")
par(xaxs="i")
image(t(remat), axes=F, col=colred)
mtext(tagline, side=3, line=.25, cex=.65)
# Add location:
rn = min(df$start)
sq = seq(-rn, rn, length.out=5)
loc = (sq - rn) / (25 * ncol(mat))
axis(side=1, at=loc, labels=FALSE)
axislabels = as.character(sapply(sq, function(x){sprintf("%0.0f",x)}))
text(x=loc, y=par()$usr[3]-0.025*(par()$usr[4] - par()$usr[3]),
     labels=axislabels, xpd=TRUE, cex=.5)
# Box around ones above 2 + add %? 
pos = par()$usr[4] - diff(par()$usr[3:4]) * frac2
abline(h=pos, lty='dotted', lwd=3)
text(x=mean(par()$usr[1:2]), y=mean(c(pos, par()$usr[4])),
     labels=label, cex=1.5, col=rgb(0,0,0,.7))
dev.off()

print("Plotting finished.")
