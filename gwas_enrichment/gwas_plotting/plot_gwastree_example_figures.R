#!/usr/bin/R
# -----------------------------------------------------------
# Script to plot/explore large + small (precomputed) examples
# NOTE: Updated 08/15/20 to plot both FDR levels into PNG
# -----------------------------------------------------------
domain = system("hostname -d", intern=TRUE)
if (length(domain) == 0) domain = ''
if (domain == 'broadinstitute.org'){
    bindir='~/data/EPIMAP_ANALYSIS/bin/'
} else {
    bindir='~/EPIMAP_ANALYSIS/bin/'
}
source(paste0(bindir, 'general_EPIMAP_ANALYSIS.R'))
library(qvalue)

# Arguments for loading data:
usetree = 'enhancers'
tol = 2500
singlematch = FALSE
plotting.only = FALSE  # No, need to run the regressions
use.adj = TRUE
use.strict = FALSE
# use.onecutoff = FALSE
use.onecutoff = TRUE
MINP=3

# Load in + process all of the relevant matrices/datasets:
commandArgs <- function(trailingOnly=TRUE){
    c(usetree, tol, singlematch, plotting.only, use.adj, use.strict, use.onecutoff) }
source(paste0(bindir, 'load_statistics_gwastree_enrichments.R'))

rm(regwide)
rm(cdll)
rm(enhdf)
gc()

# ------------
# Directories:
# ------------
gtdir = "gwas_tree_analysis/"
regdir = paste0(gtdir, "regressions/")
perdir = paste0(gtdir, "permuted_catalogs/")
epref = paste0(usetree, '_e', tol, '_')
regpref = paste0(regdir, epref)
perpref = paste0(perdir, epref)
eximgpref = sub('analysis/', 'analysis/examples/', treeimgpref)
exdir = sub('enhancers.*','', eximgpref)
cmd = paste('mkdir -p', gtdir, regdir, perdir, exdir)
system(cmd)

# Get the leaf reduced information (bag of words, max 3 terms):
leafrep = sapply(1:NN, function(x, max.terms=3){
                     blacklist = c('GLAND', 'TISSUE', 'CELL')
                     x = sub("[0-9]*_","",declist$dec[[x]])
                     # x = unique(x)
                     nx = length(x)
                     if (nx > 1){
                         term_words <- strsplit(x, "[ _,.]");
                         tab_all <- sort(table(unlist(term_words)));
                         tab_all = tab_all[tab_all > 1]
                         tab_all = tab_all[!(names(tab_all) %in% blacklist)]
                         x = paste0(tolower(head(names(sort(tab_all, decreasing=T)), max.terms)), collapse=', ')
                     } else { x = tolower(x) }
                     x = capitalize(x)
                     return(x)})
lind = which(leafrep == '')
nodetissue = nodetissue[order(nodetissue$node),]
leafrep[lind] = nodetissue$GROUP[lind]


# -------------------------------
# Plotting specific (large) GWAS:
# -------------------------------
# Choose one gwas: 
# Example of a complex trait:
strait = "Alzheimer's disease or family history of Alzheimer's disease"
sgw = aggregate(pValue ~ uid + trait, gwdf[gwdf$trait == strait,], length)

# Example of a straight-forward trait:
strait = 'Age-related macular degeneration'
suid = '23326517 - Age-related macular degeneration'
strait = "Alzheimer's disease or family history of Alzheimer's disease"
suid = "30617256 - Alzheimer's disease or family history of Alzheimer's disease"
suid = "29777097 - Alzheimer's disease or family history of Alzheimer's disease"
suid = "26192919 - Ulcerative colitis"
suid = '23326517 - Age-related macular degeneration'
suid = "26343387 - Myocardial infarction"

suid = "30617256 - Alzheimer's disease or family history of Alzheimer's disease"

suid = "26192919 - Ulcerative colitis"

# Load suffixes:
# for (suffix in c("_adj1000_1.Rda", "_adj1000_5.Rda", "_adj1000_10.Rda")){
# if (use.adj){ suffix = '_adj1000_10.Rda' } else { suffix = '.Rda' }
# if (use.strict){ suffix = '_adj1000_1.Rda' }
# Load regression for suid:
plotind = which(uids == suid)
ipref = paste0(apref, sprintf("_%05d", plotind))
print(paste(suid, '----- output with prefix', ipref))
regfile = paste0(regpref, ipref, '_lreg', suffix)
load(regfile)

# Fix with corrected:
ll$rawlp = all.regmat[suid,]
ll$df$rawlog10p = all.regmat[suid,]
ll$log10p = ll$rawlp
ll$log10p[ll$log10p > 12] = 12
ll$df$log10p = ll$log10p 

# plot(ll$df$log10p, ll$log10p)

# Setup dendrogram:
dd  = setup_dendrogram(dend3, ll, udf, declist=declist, bcutoff=MINP, altline=3)
# TODO: Update the leaf p-values. 
# TODO: Make a better labeling system for signif nodes.
# TODO: Update the legend.
# For wider/less wide (in a circle)???
cc = sqrt(cdlenlist$cons)
dd$dend = set(dd$dend, 'branches_lwd', cc / 200)

# pdf(paste0(treeimgpref, ipref, '_logreg_example_large.pdf'),width=14.5,height=12, onefile=T)
png(paste0(treeimgpref, ipref, '_logreg_example_large.png'),res=450, units='in', width=14.5,height=12)
plot.new()
circle_size = unit(1, "snpc") # snpc unit gives you a square region
pushViewport(viewport(x = 0, y = 0.5, width = circle_size, height = circle_size,
                      just = c("left", "center")))
par(omi = gridOMI(), new = TRUE)
circleplot(dend=dd$dend, lab=lab, udf=dd$udf, fractional=FALSE,
           hit.track=dd$hitdf, curved.lines=FALSE,
           leaf.pval=leafmat[strait,])
upViewport()
circos.clear()
draw(pd.legend.ext, x = circle_size, just = "left")
title(paste(suid, '-', ll$title))
dev.off()


# ---------------------------
# Plotting many (small) GWAS:
# ---------------------------
suidlist = c("26192919 - Ulcerative colitis",
             "28067908 - Crohn's disease",
             "28067908 - Inflammatory bowel disease",
             "26343387 - Myocardial infarction",
             "29212778 - Coronary artery disease",
             "30595370 - Cardiovascular disease",
             '23326517 - Age-related macular degeneration',
             "30535121 - Macular thickness",
             "30054594 - Glaucoma",
             "30617256 - Alzheimer's disease or family history of Alzheimer's disease",
             "30679814 - QT interval",
             "30275531 - Total cholesterol levels",
             "30529582 - Colorectal cancer",
             "29059683 - Breast cancer",
             "29892016 - Prostate cancer")


# List 2 :
slist2 = c('30275531 - HDL cholesterol',
           '28334899 - HDL cholesterol levels',
           '30275531 - Triglycerides',
           '28334899 - LDL cholesterol levels',
           '29769521 - Heart rate response to exercise',
           '23969696 - Fibrinogen', 
           '29422604 - Pancreatic cancer',
           '30529582 - Colorectal cancer', 
           '26831199 - Chronic kidney disease', 
           '30504769 - Gallstone disease', 
           '29403010 - Serum uric acid levels', 
           '25631615 - Optic cup area', 
           '30367059 - Thyroid stimulating hormone',
           '30586737 - Factor VIII levels', 
           '29403010 - Hemoglobin', 
           '21909110 - Blood pressure', 
           '22581228 - Fasting blood glucose')




# Run or load regression:
# # if (use.adj){ suffix = '_adj1000_10.Rda' } else { suffix = '.Rda' }
# # if (use.strict){ suffix = '_adj1000_1.Rda' }
# suffix = '_adj1000_10.Rda' 
suffix = '_adj1000_1.Rda' 
pdf(paste0(treeimgpref,apref, '_examples_small2.pdf'), width=3.5,height=4, onefile=T)
MINP=3
for (suid in slist2) {
# for (suid in suidlist) {
    # Load regression for suid:
    plotind = which(uids == suid)
    strait = unique(gwdf$trait[gwdf$uid == suid])
    ipref = paste0(apref, sprintf("_%05d", plotind))
    print(paste(suid, '----- output with prefix', ipref))
    regfile = paste0(regpref, ipref, '_lreg', suffix)
    load(regfile)
    if (sum(ll$log10p ) > 0){
        # Setup dendrogram:
        dd  = setup_dendrogram(dend3, ll, udf, declist=declist, bcutoff=MINP, altline=0)
        dd$dend = try(set(dd$dend, "branches_lwd", .5))
        # Plot small dendrogram
        NTOP=5
        ldf = ll$df
        ldf = ldf[order(ldf$log10p <= 0, ldf$pout),]
        ldf = ldf[1:NTOP,] # TOP N only
        ldf = ldf[ldf$rawlog10p > MINP,]
        ntdf = merge(ldf, nodetissue)
        ntdf = ntdf[,c('node','GROUP','COLOR')]
        names(ntdf) = c('node','symbol','color')
        ntdf = merge(ntdf, nodedf)
        ntdf$cex=.5
        par(mar=c(1,0,1,0))
        circleplot(dend=dd$dend, lab=lab, fractional=FALSE, hit.track=dd$hitdf, 
                   plot.labels=FALSE, add.metadata=TRUE, reduced.metadata=TRUE, nodedf=ntdf, scale=2)
        circos.clear()
        mtext(strait, side=3, line=0, cex=.7)
    }
}
dev.off()


# All Znam - one GWAS per trait
load(paste0(eprefix, 'kept_traits_ordered', suffix))
pdf(paste0(treeimgpref,apref, '_examples_all_small.pdf'), width=3.5,height=4, onefile=T)
MINP=3
for (suid in rev(Znam)) {
    # for (suid in suidlist) {
    # Load regression for suid:
    plotind = which(uids == suid)
    strait = unique(gwdf$trait[gwdf$uid == suid])
    ipref = paste0(apref, sprintf("_%05d", plotind))
    print(paste(suid, '----- output with prefix', ipref))
    regfile = paste0(regpref, ipref, '_lreg', suffix)
    load(regfile)
    if (sum(ll$log10p ) > 0){
        # Setup dendrogram:
        dd  = setup_dendrogram(dend3, ll, udf, declist=declist, bcutoff=MINP, altline=0)
        dd$dend = try(set(dd$dend, "branches_lwd", .5))
        # Plot small dendrogram
        NTOP=5
        ldf = ll$df
        ldf = ldf[order(ldf$log10p <= 0, ldf$pout),]
        ldf = ldf[1:NTOP,] # TOP N only
        ldf = ldf[ldf$rawlog10p > MINP,]
        ntdf = merge(ldf, nodetissue)
        ntdf = ntdf[,c('node','GROUP','COLOR')]
        names(ntdf) = c('node','symbol','color')
        ntdf = merge(ntdf, nodedf)
        ntdf$cex=.5
        par(mar=c(1,0,1,0))
        circleplot(dend=dd$dend, lab=lab, fractional=FALSE, hit.track=dd$hitdf, 
                   plot.labels=FALSE, add.metadata=TRUE, reduced.metadata=TRUE, nodedf=ntdf, scale=2)
        circos.clear()
        mtext(strait, side=3, line=0, cex=.7)
    }
}
dev.off()



# For plotting at the 0.1% FDR level:
pcutdf = read.delim(paste0(gtdir, 'cutoffs_l10', setpref, '_all_pernode.tsv'))
fper1p = pcutdf[pcutdf$type == 'All 0.1%','cut']
rmatper1p = all.regmat / fper1p
all.regmat1p = 1 * (rmatper1p > 1) * all.regmat
sum(apply(all.regmat1p[guid,], 1, max) > 0)


# -----------------------------------------------------
# All Znam in all GWAS (used for construction of plots)
# -----------------------------------------------------
# NOTE: Missing before: 26920376 - Multiple sclerosis and body mass index (pleiotropy)
load(paste0(extpref, 'kept_allgwas_ordered', cutsuff, suffix))
print(length(Znam))
# Write out the 0.1% list:
Zfix = Znam 
Zfix[Zfix == "26974007 - Chronic inflammatory diseases (pleiotropy)"] = paste0("26974007 - Chronic inflammatory diseases (ankylosing spondylitis, ",
                                                                       "Crohn's disease, psoriasis, primary sclerosing cholangitis, ulcerative colitis) (pleiotropy)")
Znam1p = Zfix[which(apply(all.regmat1p[Zfix,], 1, max) > 0)]
write(rev(Znam1p), 'tree_gwas1p_list.txt')

MINP=3
for (j in 1:length(Znam)){
    suid = rev(Znam)[j]
    # Load regression for suid:
    if (suid == "26974007 - Chronic inflammatory diseases (pleiotropy)"){
        suid = paste0("26974007 - Chronic inflammatory diseases (ankylosing spondylitis, ",
                      "Crohn's disease, psoriasis, primary sclerosing cholangitis, ulcerative colitis) (pleiotropy)")
    } 
    plotind = which(uids == suid)
    strait = unique(gwdf$trait[gwdf$uid == suid])
    spmid = unique(gwdf$pubMedID[gwdf$uid == suid])
    sgwssdf = gwssdf[gwssdf$uid == suid,]
    sinit = sgwssdf$initSample[order(sgwssdf$sampsize, decreasing=T)][1]
    s2init = split.text(sinit, width=90)
    s2trait = split.text(strait, width=50)
    ipref = paste0(apref, sprintf("_%05d", plotind))
    cat(j, '\t', suid, '\t', ipref, '\t')
    regfile = paste0(regpref, ipref, '_lreg', suffix)
    # Load file, but use the adjusted vs. one-cut cutoff
    load(regfile)
    trait = sub(".* - ", "", suid)
    pmid = sub(" - .*", "", suid)
    traitstr = gsub("'","_", gsub(" ", "_", tolower(trait)))
    traitstrnoparen = gsub("/", "_", gsub("\\)" ,"", gsub("\\(","",traitstr)))

    for (pcut in c(1, 0.1)){
        # Fix with corrected at either threshold:
        if (pcut == 1){
            ll$rawlp = all.regmat[suid,]
            ll$df$rawlog10p = all.regmat[suid,]
            pcutstr = '_1pct'
        } else {
            ll$rawlp = all.regmat1p[suid,]
            ll$df$rawlog10p = all.regmat1p[suid,]
            pcutstr = '_pt1pct'
        }
        ll$log10p = ll$rawlp
        ll$log10p[ll$log10p > 12] = 12
        ll$df$log10p = ll$log10p 

        if (sum(ll$log10p) > 0){
            cat(pcut, '\t')
            # Setup dendrogram:
            dd  = setup_dendrogram(dend3, ll, udf, declist=declist, bcutoff=MINP, altline=0)
            dd$dend = try(set(dd$dend, "branches_lwd", .5))
            # Plot small dendrogram
            NTOP=10
            ldf = ll$df
            ldf = ldf[order(ldf$log10p <= 0, ldf$pout),]
            ldf = ldf[1:NTOP,] # TOP N only
            ldf = ldf[ldf$rawlog10p > MINP,]
            ldf$rank = 1:nrow(ldf)
            ntdf = merge(ldf, nodetissue)
            ntdf$GROUP = paste0(ntdf$rank,'. ', ntdf$GROUP, "\n(", leafrep[ntdf$node], ")")
            ntdf = ntdf[,c('node','GROUP','COLOR')]
            names(ntdf) = c('node','symbol','color')
            ntdf = merge(ntdf, nodedf)
            ntdf$cex=.4

            imgfilename = paste0(eximgpref, apref, '_gwas_small_', traitstrnoparen, '_', pmid, pcutstr, '.png')
            if (!file.exists(imgfilename)){
                png(imgfilename, res=400, units='in', width=3.5,height=4)
                par(mar=c(1,0,1,0))
                circleplot(dend=dd$dend, lab=lab, fractional=FALSE, hit.track=dd$hitdf, 
                           plot.labels=FALSE, add.metadata=TRUE, reduced.metadata=TRUE, nodedf=ntdf, scale=2)
                circos.clear()
                if (s2trait != strait){
                    mtext(s2trait, side=3, line=-.5, cex=.5)
                    mtext(spmid, side=3, line=-1, cex=.3)
                } else {
                    mtext(s2trait, side=3, line=0, cex=.7)
                    mtext(spmid, side=3, line=-.75, cex=.5)
                }
                mtext(s2init, side=1, line=0, cex=.3)
                dev.off()
            }
        }
    }
    cat('Done\n')

}
# dev.off()


# ----------------------------------------------
# All Znam in all GWAS past certain sample size:
# ----------------------------------------------
load(paste0(eprefix, 'kept_allgwas_ordered', suffix))
pdf(paste0(treeimgpref,apref, '_examples_all_gwas_small_50k.pdf'), width=3.5,height=4, onefile=T)
MINP=3
for (suid in rev(Znam)) {
    # for (suid in suidlist) {
    # Load regression for suid:
    plotind = which(uids == suid)
    strait = unique(gwdf$trait[gwdf$uid == suid])
    spmid = unique(gwdf$pubMedID[gwdf$uid == suid])
    sgwssdf = gwssdf[gwssdf$uid == suid,]
    ssamp = sgwssdf$sampsize[order(sgwssdf$sampsize, decreasing=T)][1]
    if (ssamp > 50000){
        sinit = sgwssdf$initSample[order(sgwssdf$sampsize, decreasing=T)][1]
        s2init = split.text(sinit, width=90)
        ipref = paste0(apref, sprintf("_%05d", plotind))
        print(paste(suid, '----- output with prefix', ipref))
        regfile = paste0(regpref, ipref, '_lreg', suffix)
        load(regfile)
        if (sum(ll$log10p) > 0){
            # Setup dendrogram:
            dd  = setup_dendrogram(dend3, ll, udf, declist=declist, bcutoff=MINP, altline=0)
            dd$dend = try(set(dd$dend, "branches_lwd", .5))
            # Label the top 5 nodes with tissue: 
            NTOP=5
            ldf = ll$df
            ldf = ldf[order(ldf$log10p <= 0, ldf$pout),]
            ldf = ldf[1:NTOP,] # TOP N only
            ldf = ldf[ldf$rawlog10p > MINP,]
            ntdf = merge(ldf, nodetissue)
            ntdf = ntdf[,c('node','GROUP','COLOR')]
            names(ntdf) = c('node','symbol','color')
            ntdf = merge(ntdf, nodedf)
            ntdf$cex=.5
            # Plot small dendrogram
            par(mar=c(1,0,1,0))
            circleplot(dend=dd$dend, lab=lab, fractional=FALSE, hit.track=dd$hitdf, 
                       plot.labels=FALSE, add.metadata=TRUE, reduced.metadata=TRUE, nodedf=ntdf, scale=2)
            circos.clear()
            mtext(strait, side=3, line=0, cex=.7)
            mtext(spmid, side=3, line=-.75, cex=.5)
            mtext(s2init, side=1, line=0, cex=.3)
        }
    }
}
dev.off()



# -----------------------------------
# Plot dendrogram - with genes on it!
# -----------------------------------
# Temporary:
# labels(dd$dend) = paste0(labels(dd$dend),'_',labels(dend))
pdf(paste0(treeimgpref,apref, '_logreg_test.pdf'),width=14.5,height=12, onefile=T)
plot.new()
circle_size = unit(1, "snpc") # snpc unit gives you a square region
pushViewport(viewport(x = 0, y = 0.5, width = circle_size, height = circle_size,
                      just = c("left", "center")))
par(omi = gridOMI(), new = TRUE)
circleplot(dend=dd$dend, lab=lab, udf=dd$udf, fractional=FALSE, hit.track=dd$hitdf, leaf.pval=leafmat[strait,])
upViewport()
circos.clear()
# draw(pd.legend.ext, x = circle_size, just = "left")
title(paste(suid, '-', ll$title))
dev.off()



