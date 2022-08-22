#!/usr/bin/R
# ---------------------------
# Plot + rank the IC metrics:
# ---------------------------
domain = system("hostname -d", intern=TRUE)
if (length(domain) == 0) domain = ''
if (domain == 'broadinstitute.org'){
    bindir='~/data/EPIMAP_ANALYSIS/bin/'
} else {
    bindir='~/EPIMAP_ANALYSIS/bin/'
}
source(paste0(bindir, 'general_EPIMAP_ANALYSIS.R'))

library(ggplot2)
library(ggpubr)
library(gridExtra)

# ------------------------------
# Load metrics and mapping data:
# ------------------------------
method = '_qnorm'
# method = '_simple'
df = read.delim(paste0('aggregated_IC_stats',method,'.tsv'), header=T, sep="\t")
df = gather(df, metric, value, -team, -track)
df$ic.id = sub("M.*","", df$track)
df$ic.mark = sub(".*M","M", df$track)

nmap = read.delim('Annotation/ic_namespace_mapping.tsv',header=T)
mmap = read.delim('Annotation/ic_mark_mapping.tsv',header=F)
names(mmap) = c('ic.mark', 'mark')
tmap = read.delim('Annotation/team_name_round2.tsv', header=F)
names(tmap) = c('team', 'teamname')
tmap$teamname = sub("Hongyang Li and Yuanfang Guan","HLYG", tmap$teamname)
df = merge(merge(df, mmap, all.x=TRUE), nmap, all.x=TRUE)


# ------------------------------------------------------------------
# Score the teams according to the IC procedure, without bootstraps:
# ------------------------------------------------------------------
df$rank = NA
df = df[order(df$track),]
for (track in unique(df$track)){
    for (metric in unique(df$metric)){
        ind = which(df$track == track & df$metric == metric)
        score = df$value[ind]
        rank = order(order(-score))
        df$rank[ind] = rank
    }
}

# Re-rank tracks:
aggdf = aggregate(rank ~ track + team, df, mean)
aggdf = aggdf[order(aggdf$track),]
aggdf$track.rank = NA
for (track in unique(df$track)){
    ind = which(aggdf$track == track)
    score = aggdf$rank[ind]
    aggdf$track.rank[ind] = order(order(score))
}

# Aggregate the track ranks:
# NOTE: Not using the min(0.5, r_e), unclear how that works if rankings are 1-27.
rankdf = aggregate(track.rank ~ team, aggdf, mean) 
rankdf$final.rank = order(order(rankdf$track.rank))
rankdf = merge(rankdf, tmap)
rankdf = rankdf[order(rankdf$final.rank),]
rownames(rankdf) = NULL

png(paste0("IC_rankings_relativemetrics",method,".png"), res=250, units='in', width=6, height=8)
grid.table(rankdf)
dev.off()

pdf(paste0("IC_rankings_relativemetrics",method,".pdf"), width=6, height=8)
grid.table(rankdf)
dev.off()


# Plot the track.ranks on each experiment:
# ----------------------------------------
aggdf$ic.id = sub("M.*","", aggdf$track)
aggdf$ic.mark = sub(".*M","M", aggdf$track)
aggdf = merge(merge(aggdf, mmap, all.x=TRUE), nmap, all.x=TRUE)
aggdf = merge(aggdf, tmap)
aggdf$teamname = factor(aggdf$teamname, levels=rankdf$teamname)

colvals = rep(c('black','grey75'), 13)
names(colvals) = rankdf$teamname
colvals['Avocado'] = 'forestgreen'
colvals['Average'] = 'orange'
colvals['Lavawizard'] = 'red'
colvals['imp'] = 'blue'
colvals['HLYG v1'] = 'purple'
colvals['ChromImpute'] = 'goldenrod1'

gp = ggplot(aggdf, aes(teamname, track.rank, color=factor(teamname))) + 
    facet_wrap(~mark, scales='free_y') + 
    theme_pubr() + 
    scale_color_manual(values=colvals) +
    labs(x='Team', y='Track Rank') + 
    geom_jitter(cex=.8) + 
    theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5)) + 
    theme(legend.position='none')
ggsave(paste0('IC_rankings_trackranks', method,'.png'), dpi=350, units='in', width=13, height=8)
ggsave(paste0('IC_rankings_trackranks', method,'.pdf'), dpi=350, units='in', width=13, height=8)

resdf = merge(aggdf, rankdf[,c('teamname','final.rank')])
resdf$rank.residual = resdf$track.rank - resdf$final.rank
gp = ggplot(resdf, aes(teamname, rank.residual, color=factor(teamname))) + 
    facet_wrap(~mark, scales='free_y') + 
    theme_pubr() + 
    scale_color_manual(values=colvals) +
    geom_hline(yintercept=0, lty='dashed') + 
    geom_jitter(cex=.8) + 
    labs(x='Team', y='Track Rank - Team Rank') + 
    theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5)) + 
    theme(legend.position='none')
ggsave(paste0('IC_rankings_trackrank_residuals',method,'.png'), dpi=350, units='in', width=13, height=8)
ggsave(paste0('IC_rankings_trackrank_residuals',method,'.pdf'), dpi=350, units='in', width=13, height=8)


# ----------------------
# Plot the actual values
# ----------------------
colvals2 = colvals
colvals2[colvals2 == 'black'] = 'grey75'
df = merge(df, tmap)
ggplot(df, aes(mark, value, color=factor(teamname))) + 
    facet_wrap(~metric, scales='free_y') + 
    scale_color_manual(values=colvals2) + 
    theme_pubr() + 
    geom_boxplot(outlier.shape=NA, color='black') +
    geom_jitter()


