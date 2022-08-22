#!/usr/bin/R
# ----------------------------------------------------------
# Compute imputed agreement metrics for chr1 on the IC data:
# ----------------------------------------------------------
domain = system("hostname -d", intern=TRUE)
if (length(domain) == 0) domain = ''
if (domain == 'broadinstitute.org'){
    bindir='~/data/EPIMAP_ANALYSIS/bin/'
} else {
    bindir='~/EPIMAP_ANALYSIS/bin/'
}
source(paste0(bindir, 'general_EPIMAP_ANALYSIS.R'))

library(PRROC)

# Read in prefix:
args = (commandArgs(TRUE))
if (length(args) == 0) {
    print("No arguments given")
} else {        
    prefix = args[1]
    infofile = args[2]
}

baseprefix = sub(".*/", "", prefix)
calc.25 = FALSE

print('[STATUS] Using the following prefix + base:')
print(prefix)
print(baseprefix)

# --------------------------------------
# Read in datasets imputed and external:
# --------------------------------------
info = read.delim(infofile, header=F)
names(info) = c('id','mark', 'fpref')

# Read in 25bp and 200bp data:
df200 = read.delim(gzfile(paste0(prefix, '.chr1.b200.lo.final.bed.gz')), header=F, stringsAsFactors=F)
names(df200) = c('chr','start','end', 'value')
ind200 = which(df200$value != 0) # Do not compare on mis-mapped areas.
or200 = df200$value[ind200]
ob200 = 1 * (or200 >= 2)

# 25bp data:
df25 = read.delim(gzfile(paste0(prefix, '.chr1.b25.lo.final.bed.gz')), header=F, stringsAsFactors=F)
names(df25) = c('chr','start','end', 'value')
ind25 = which(df25$value != 0) # Do not compare on mis-mapped areas.
or25 = df25$value[ind25]
ob25 = 1 * (or25 >= 2)


BinMean <- function (vec, every, na.rm = FALSE) {
    n <- length(vec)
    x <- .colMeans(vec, every, n %/% every, na.rm)
    r <- n %% every
    if (r) x <- c(x, mean.default(vec[(n - r + 1):n], na.rm = na.rm))
    x
}

comp.io = function(or, ob, ir){
    # Sort scores for faster roc, etc:
    ord = order(ir)
    io = ir[ord]
    ro = or[ord]
    oo = ob[ord]
    # Perform correlation + prediction comparisons:
    gwcorr = cor(or, ir)
    cat('CORR:', '\t', round(gwcorr,3), '\n')
    roc = roc.curve(scores.class0=io[oo == 1],
                    scores.class1=io[oo == 0], curve=TRUE, sorted=TRUE)
    pr = pr.curve(scores.class0=io[oo == 1],
                  scores.class1=io[oo == 0], curve=TRUE, sorted=TRUE)
    cat('AUC:', '\t', round(roc$auc,3), '\t', round(pr$auc.integral,3), '\n')
    # ROC and AUC on the top 1%:
    thresh = quantile(ro, .99)
    oo2 = 1 * (ro >= thresh)
    roc2 = roc.curve(scores.class0=io[oo2 == 1],
                    scores.class1=io[oo2 == 0], curve=TRUE, sorted=TRUE)
    pr2 = pr.curve(scores.class0=io[oo2 == 1],
                  scores.class1=io[oo2 == 0], curve=TRUE, sorted=TRUE)
    cat('AUC1:', '\t', round(roc2$auc,3), '\t', round(pr2$auc.integral,3), '\n')
    # Calculate the catch metrics:
    ot1 = 1 * (or >= quantile(or, .99))
    ot5 = 1 * (or >= quantile(or, .95))
    it1 = 1 * (ir >= quantile(ir, .99))
    it5 = 1 * (ir >= quantile(ir, .95))
    # Catch top 1 % from top 5%
    co1 = sum(ot1 * it5) / sum(ot1)
    ci1 = sum(it1 * ot5) / sum(it1)
    mi1 = sum(it1 * ot1) / sum(it1)
    mo1 = sum(it1 * ot1) / sum(ot1)
    # Aggregate the statistics:
    sdf = data.frame(gwcorr=gwcorr,
                     auroc=roc$auc, auprc=pr$auc.integral,
                     auroc1=roc2$auc, auprc1=pr2$auc.integral,
                     catchobs=co1, catchimp=ci1,
                     matchobs=mo1, matchimp=mi1)
    return(sdf)
}


# ------------------------------------------------------------
# Read imputed track in compare to the external observed data:
# ------------------------------------------------------------
x25 = as.numeric(scan(paste0(prefix, '.chr1.imputed.b25.raw.txt.gz'),'c'))
x200 = BinMean(x25, 8)
x25 = x25[ind25]
x200 = x200[ind200]

# Compute their statistics:
statsdf = comp.io(or200, ob200, x200)
statsdf$step = 200
statsdf$comp.id = baseprefix

if (calc.25){
    sdf = comp.io(or25, ob25, x25)
    sdf$step = 25
    sdf$comp.id = baseprefix
    statsdf = rbind(statsdf, sdf)
}

# -------------------------------------------------
# Read in the observed data and compare one by one:
# -------------------------------------------------
m25 = NULL
print(nrow(info))
for (i in 1:nrow(info)){
    id = as.character(info$id[i])
    print(id)
    fname = paste0('ChromImpute/imputed/chr1_', info$fpref[i], '.wig.gz')
    x25 = scan(fname,'c')
    x25 = as.numeric(x25[3:length(x25)])

    # Add to mean counter:
    if (is.null(m25)){
        m25 = x25
        j = 1
    } else {
        m25 = m25 + x25
        j = j + 1 
    }
    
    # Compute the 200bp mean and reduce to comparable indices:
    x200 = BinMean(x25, 8)
    x25 = x25[ind25]
    x200 = x200[ind200]

    # Compute metrics:
    sdf = comp.io(or200, ob200, x200)
    sdf$step = 200
    sdf$comp.id = id 
    statsdf = rbind(statsdf, sdf)

    # 25bp
    if (calc.25){
        sdf = comp.io(or25, ob25, x25)
        sdf$step = 25
        sdf$comp.id = id
        statsdf = rbind(statsdf, sdf)
    }
}

# Compute for the mean as well:
x25 = m25 / j
x200 = BinMean(x25, 8)
x25 = x25[ind25]
x200 = x200[ind200]

# Compute metrics:
sdf = comp.io(or200, ob200, x200)
sdf$step = 200
sdf$comp.id = 'Mean'
statsdf = rbind(statsdf, sdf)

# 25bp
if (calc.25){
    sdf = comp.io(or25, ob25, x25)
    sdf$step = 25
    sdf$comp.id = 'Mean'
    statsdf = rbind(statsdf, sdf)
}

# Add prefix, write out:
statsdf$prefix = baseprefix
write.table(statsdf, file=paste0(prefix, '.chr1.stats.tsv'), quote=F, row.names=F, sep="\t")

