#!/usr/bin/R
# ---------------------------------------------------
# NOTE: original version can be found at count_snps_tree.R
# Preliminary analysis for tree GWAS - 
# Counts on a fixed matrix from epigenomic similarity
# ---------------------------------------------------
domain = system("hostname -d", intern=TRUE)
if (length(domain) == 0) domain = ''
if (domain == 'broadinstitute.org'){
    bindir='~/data/EPIMAP_ANALYSIS/bin/'
} else {
    bindir='~/EPIMAP_ANALYSIS/bin/'
}
source(paste0(bindir, 'general_EPIMAP_ANALYSIS.R'))
source(paste0(bindir, 'auxiliary_chromImpute_functions.R'))

# Arguments:
usetree = 'enhancers'
tol = 2500
singlematch = FALSE # Use 1-to-1 SNP enhancer mapping only?
plotting.only = FALSE

# Load in + process all of the relevant matrices/datasets:
commandArgs <- function(trailingOnly=TRUE){
    c(usetree, tol, singlematch, plotting.only) }
source(paste0(bindir, 'load_gwastree_analysis.R'))

print(paste("Images go to: ", treeimgpref))

# ----------------
# Plot basic tree:
# ----------------
pdf(paste0(treeimgpref, 'basic_tree.pdf'), width=14.5, height=12)
plot.new()
circle_size = unit(1, "snpc") # snpc unit gives you a square region
pushViewport(viewport(x = 0, y = 0.5, width = circle_size, height = circle_size,
    just = c("left", "center")))
par(omi = gridOMI(), new = TRUE)
circleplot(dend3, lab=lab, udf=udf)
upViewport()
draw(pd.legend, x = circle_size, just = "left")
dev.off()

# ----------------------
# Look at specific GWAS:
# ----------------------
# List of traits:
ut = unique(gwdf$trait)
# Choose one gwas: 
# strait = "Crohn's disease"
strait = 'Glomerular filtration rate'
strait = 'LDL cholesterol'
strait = 'Macular thickness'
strait = 'QT interval'

ll = get_pvalues(strait, dflist, cdlenlist=cdlenlist, pdf=pdf,
                 ingroup='cons', outgroup='enh', against.parent=FALSE)
# ldf = ll$df[order(ll$df$pout),]
# head(ldf,50)
dd  = setup_dendrogram(dend3, ll, udf, declist=declist)

# Plot basic dendrogram
pdf(paste0(treeimgpref,'full_test.pdf'),width=14.5,height=12, onefile=T)
plot.new()
circle_size = unit(1, "snpc") # snpc unit gives you a square region
pushViewport(viewport(x = 0, y = 0.5, width = circle_size, height = circle_size,
                      just = c("left", "center")))
par(omi = gridOMI(), new = TRUE)
circleplot(dend=dd$dend, lab=lab, udf=dd$udf, fractional=FALSE)
upViewport()
circos.clear()
draw(pd.legend.ext, x = circle_size, just = "left")
title(paste(strait, '-', ll$title))
dev.off()


# -----------------------------
# Logistic regression approach:
# -----------------------------

# For testing, choose relatively easy trait:
# strait = 'LDL cholesterol'
strait = 'Macular thickness'
testtraits = c('Macular thickness', "QT interval", 'LDL cholesterol', "Crohn's disease")
strait = 'LDL cholesterol'
strait = "Crohn's disease"
strait = "QT interval"
type = 'cons'  # Also do union and average
against = 'parent'  # Also do sibling and all enhancers
weighted = FALSE
# weighted = TRUE  # Doesn't work currently. Check math.

# Process all traits x types x against 
for (type in c('cons','union')){
    against = 'parent'  # Also do sibling and all enhancers
    apref = paste0(type, '_', against)
    if (weighted){ 
        weights = sqrt(1 / matmarg[,2])
        apref = paste0(apref, '_weighted')
    } else {
        weights = NULL
    }
    etfile = paste0('logreg_', apref, midpref, 'pvals_elist.Rda')
    if (!file.exists(etfile)){
        lrlist = list()
        for (strait in expandedlist){
            ll = get_disjoint_lr(strait, type=type, against=against, weights=weights)
            print(paste(strait, ' -- ', sum(ll$log10p > 0)))
            lrlist[[strait]] = ll
        }
        save(lrlist, file=etfile)
    } else {
        load(etfile)
    }


    # -------------------------------------------------------
    # Plot basic dendrogram - with more highlighted branches:
    # -------------------------------------------------------
    # strait = 'Systolic blood pressure'
    ll = lrlist[[strait]]
    dd  = setup_dendrogram(dend3, ll, udf, declist=declist, bcutoff=3)
    pdf(paste0(treeimgpref,apref, '_logreg_test.pdf'),width=14.5,height=12, onefile=T)
    plot.new()
    circle_size = unit(1, "snpc") # snpc unit gives you a square region
    pushViewport(viewport(x = 0, y = 0.5, width = circle_size, height = circle_size,
                          just = c("left", "center")))
    par(omi = gridOMI(), new = TRUE)
    circleplot(dend=dd$dend, lab=lab, udf=dd$udf, fractional=FALSE, hit.track=dd$hitdf, leaf.pval=leafmat[strait,])
    upViewport()
    circos.clear()
    draw(pd.legend.ext, x = circle_size, just = "left")
    title(paste(strait, '-', ll$title))
    dev.off()

    # -----------------------------
    # Plot all expanded trait list:
    # -----------------------------
    pdf(paste0(treeimgpref, apref, '_logreg_examples_bcut3.pdf'),width=14.5,height=12, onefile=T)
    for (strait in expandedlist){
        print(strait)
        ll = lrlist[[strait]]
        dd  = setup_dendrogram(dend3, ll, udf, declist=declist, bcutoff=3)
        plot.new()
        circle_size = unit(1, "snpc") # snpc unit gives you a square region
        pushViewport(viewport(x = 0, y = 0.5, width = circle_size, height = circle_size, just = c("left", "center")))
        par(omi = gridOMI(), new = TRUE)
        circleplot(dend=dd$dend, lab=lab, udf=dd$udf, 
                   fractional=FALSE, add.metadata=FALSE, hit.track=dd$hitdf, leaf.pval=leafmat[strait,] )
        upViewport()
        circos.clear()
        # draw(pd.legend.ext, x = circle_size, just = "left")
        title(paste(strait, '-', ll$title))
    }
    dev.off()

    pdf(paste0(treeimgpref, apref, '_logreg_examples_bcut6.pdf'),width=14.5,height=12, onefile=T)
    for (strait in expandedlist){
        print(strait)
        ll = lrlist[[strait]]
        dd  = setup_dendrogram(dend3, ll, udf, declist=declist, bcutoff=6)
        plot.new()
        circle_size = unit(1, "snpc") # snpc unit gives you a square region
        pushViewport(viewport(x = 0, y = 0.5, width = circle_size, height = circle_size, just = c("left", "center")))
        par(omi = gridOMI(), new = TRUE)
        circleplot(dend=dd$dend, lab=lab, udf=dd$udf, fractional=FALSE,
                   add.metadata=FALSE, hit.track=dd$hitdf, leaf.pval=leafmat[strait,])
        upViewport()
        circos.clear()
        # draw(pd.legend.ext, x = circle_size, just = "left")
        title(paste(strait, '-', ll$title))
    }
    dev.off()
}

# ------------------------------------------
# Compute and plot all traits for consensus:
# ------------------------------------------
alltraits = unique(gwdf$trait)
all.tfile = paste0('logreg_', apref, midpref, 'pvals_alltrait.Rda')
type ='cons'
against = 'parent'  # Also do sibling and all enhancers
apref = paste0(type, '_', against)
if (!file.exists(all.tfile)){
    lrlist = list()
    for (strait in alltraits){
        ll = get_disjoint_lr(strait, type=type, against=against)
        print(paste(strait, ' -- ', sum(ll$log10p > 0)))
        lrlist[[strait]] = ll
    }
    save(lrlist, file=all.tfile)
} else {
    load(all.tfile)
}


ut = unique(gwdf$trait)
# NOTE: ONLY PLOT IF SIGNIF?





# NOTE: Fast method: Eval the model node + 1 1665 times, each time start with the parent's coefficients. Need to make sure can use the chisq statistics for this.

# - here, we have to do 1.5 regression per node (parent + node!)







# ==================================
# Plotting different sets of traits:
# ==================================
testlist = list(c('union','union',TRUE),
                # c('novel','union',TRUE),
                c('union','enh',FALSE),
                # c('novel','enh',FALSE),
                # c('diff','enh', FALSE),
                c('diff','cons',FALSE),
                c('cons','enh', FALSE))

# Pick four semi-obvious GWAS:
testtraits = c('Macular thickness', "QT interval", 'LDL cholesterol', "Crohn's disease")
pdf(paste0(treeimgpref,'full_test_types.pdf'),width=13,height=12, onefile=T)
for (strait in testtraits) {
    for (testset in testlist){
        print(testset)
        ll = get_pvalues(strait, dflist, cdlenlist=cdlenlist, pdf=pdf,
                         ingroup=testset[1], outgroup=testset[2], against.parent=testset[3])
        dd  = setup_dendrogram(dend3, ll, udf, declist=declist, altline=0)
        plot.new()
        circle_size = unit(1, "snpc") # snpc unit gives you a square region
        pushViewport(viewport(x = 0, y = 0.5, width = circle_size, height = circle_size,
                              just = c("left", "center")))
        par(omi = gridOMI(), new = TRUE)
        circleplot(dend=dd$dend, lab=lab, udf=dd$udf, add.metadata=FALSE)
        upViewport()
        circos.clear()
        # draw(pd.legend.ext, x = circle_size, just = "left")
        title(paste(strait, '-', ll$title))
    }
}
dev.off()


# ---------------------------
# Plot many different traits:
# ---------------------------

# EVAL MANY, reduced size:
pdf(paste0(treeimgpref,'full_examples_gwaslist.pdf'),width=13,height=12, onefile=T)
for (strait in tlist) {
    testset = c('cons','enh', FALSE)
    print(strait)
    ll = get_pvalues(strait, dflist, cdlenlist=cdlenlist, pdf=pdf,
                     ingroup=testset[1], outgroup=testset[2], against.parent=testset[3])
    dd  = setup_dendrogram(dend3, ll, udf, declist=declist, altline=0)  # Plot with none-lines
    plot.new()
    circle_size = unit(1, "snpc") # snpc unit gives you a square region
    pushViewport(viewport(x = 0, y = 0.5, width = circle_size, height = circle_size,
                          just = c("left", "center")))
    par(omi = gridOMI(), new = TRUE)
    circleplot(dend=dd$dend, lab=lab, udf=dd$udf, add.metadata=FALSE)  # Plot without metadata!
    upViewport()
    circos.clear()
    title(paste(strait, '-', ll$title))
}
dev.off()


NCLUST=10
# ONLY TEST diff/enh
pdf(paste0(treeimgpref,"split_k",NCLUST,"_examples.pdf"), width=14.5, height=12, onefile=T)
for (strait in expandedlist){
    print(strait)
    testlist = list(c('diff','enh', FALSE))
    for (testset in testlist){
        ll = get_pvalues(strait, dflist, cdlenlist=cdlenlist, pdf=pdf,
                         ingroup=testset[1], outgroup=testset[2], against.parent=testset[3])
        dd  = setup_dendrogram(dend3, ll, udf, declist=declist)
        plot.new()
        circle_size = unit(1, "snpc") # snpc unit gives you a square region
        pushViewport(viewport(x = 0, y = 0.5, width = circle_size, height = circle_size,
                              just = c("left", "center")))
        par(omi = gridOMI(), new = TRUE)
        circleplot(dend=dd$dend, lab=lab, udf=dd$udf, fulltrack=FALSE)
        upViewport()
        circos.clear()
        # draw(pd.legend.ext, x = circle_size, just = "left")
        title(paste(strait, '-', ll$title))
    }
}
dev.off()


# Test methods:
testlist = list(c('union','union',TRUE), 
                c('diff','enh', FALSE))
# c('novel','union',TRUE)
# c('diff','cons',FALSE)
# c('cons','enh', FALSE)
# Chunk traits:
chunksize=20
nchunk = ceiling(length(expandedlist) / chunksize)
for (i in 1:nchunk){
    pdf(paste0(treeimgpref,"examples_comp_methods", sprintf("%03d",i), ".pdf"), width=14.5, height=12, onefile=T)
    chunklist = expandedlist[((i-1) * chunksize + 1):min(i * chunksize, length(expandedlist))]
    print(i)
    print(chunklist)
    for (strait in chunklist){
        print(strait)
        for (testset in testlist){
            ll = get_pvalues(strait, dflist, cdlenlist=cdlenlist, pdf=pdf,
                             ingroup=testset[1], outgroup=testset[2], against.parent=testset[3])
            dd  = setup_dendrogram(dend3, ll, udf, declist=declist)
            plot.new()
            circle_size = unit(1, "snpc") # snpc unit gives you a square region
            pushViewport(viewport(x = 0, y = 0.5, width = circle_size, height = circle_size,
                                  just = c("left", "center")))
            par(omi = gridOMI(), new = TRUE)
            circleplot(dend=dd$dend, lab=lab, udf=dd$udf, fractional=FALSE)
            upViewport()
            circos.clear()
            draw(pd.legend.ext, x = circle_size, just = "left")
            title(paste(strait, '-', ll$title))
        }
    }
    dev.off()
}



# Look at Triglycerides - do we need to prune, and how:
strait = 'Triglycerides'
sgw = aggregate(pValue ~ uid + trait, gwdf[gwdf$trait == strait,], length)
suid = sgw[order(sgw$pValue, decreasing=T), 'uid'][1]
subdf = gwdf[gwdf$uid == suid,]
dim(subdf)
subdf = aggregate(pValue ~ chrom + chromStart + chromEnd + uid + trait, subdf, min)
dim(subdf)

subdf = subdf[order(subdf$chrom, subdf$chromStart),]

# Want to prune - take top snp in +/-1Mb at all times

pdf(paste0(treeimgpref,'full_test_types.pdf'),width=13,height=12, onefile=T)
for (testset in testlist){
    print(testset)
    ll = get_pvalues(strait, dflist, cdlenlist=cdlenlist, pdf=pdf,
                     ingroup=testset[1], outgroup=testset[2], against.parent=testset[3])
    dd  = setup_dendrogram(dend3, ll, udf, declist=declist)
    plot.new()
    circle_size = unit(1, "snpc") # snpc unit gives you a square region
    pushViewport(viewport(x = 0, y = 0.5, width = circle_size, height = circle_size,
                          just = c("left", "center")))
    par(omi = gridOMI(), new = TRUE)
    circleplot(dend=dd$dend, lab=lab, udf=dd$udf)
    upViewport()
    circos.clear()
    # draw(pd.legend.ext, x = circle_size, just = "left")
    title(paste(strait, '-', ll$title))
}
dev.off()


