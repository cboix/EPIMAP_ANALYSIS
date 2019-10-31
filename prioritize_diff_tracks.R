#!/usr/bin/R
# -----------------------------
# Prioritize impobs diff tracks
# -----------------------------
domain = system("hostname -d", intern=TRUE)
if (domain == 'broadinstitute.org'){
    source('~/data/EPIMAP_ANALYSIS/bin/general_EPIMAP_ANALYSIS.R')
} else {
    source('~/EPIMAP_ANALYSIS/bin/general_EPIMAP_ANALYSIS.R')
}
options(scipen=45) # so we dont get issues writing integers into bedfiles
today <- format(Sys.time(), "%m%d%y")
library(dplyr)
library(ggplot2)
# library(ggrepel)

# Arguments:
args=(commandArgs(TRUE))
if (length(args)==0) {
    stop("No arguments supplied: Need table filename.")
} else {        
    tsvfile = args[1]
    if (length(args) > 1){ mark = args[2] }
    if (length(args) > 2){ outfile = args[3] }
}

# Read in diffdistance data:
df = read.delim(tsvfile, header=T, stringsAsFactors=F)
df = df[df$mark == mark & df$against !=mark,]

# Predict the distance:
df$dist = df$dist / max(df$dist)
model = glm(dist ~ sample + against, df, family='gaussian')
pred = predict(model, df)

# Look for points with lower than expected:
df$diff = (df$dist - pred) / pred
df$pred = pred
df = df[order(df$diff),]
labdf=df[df$dist < 0.5,]
 
# Get mark specific quantile:
df$quant=0
for (ag in unique(df$against)){
    idx = which(df$against == ag)
    df$quant[idx] = rank(df$dist[idx])/ length(idx)
}

# Plot metric:
labdf=df[df$quant < 0.05,]
gplot = ggplot(df, aes(dist, dist-pred)) +
    geom_point(alpha=0.5, color='black') +
    facet_wrap(~against) + 
    lims(x=c(0,1), y=c(-1,1)) + 
    labs(x='Distance', y='Distance - Predicted Distance') + 
    geom_hline(yintercept=0) + 
    theme_bw() +
    geom_text(data=labdf, aes(x=dist, y=dist-pred, label=sample), size=2,vjust='outward', adj=1.1)
    # geom_text_repel(data=labdf, aes(x=dist, y=dist-pred, label=sample), size=3)
ggsave(paste0(outfile,'_metrics.png'), gplot,dpi=300,units='in',width=8,height=7)

# Write 
outdf = df[, c('sample','mark','against','diff', 'quant')]
write.table(outdf, paste0(outfile, '.tsv'), row.names=F, col.names=F, quote=F, sep="\t")
