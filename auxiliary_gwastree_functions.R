#!/usr/bin/R
# --------------------------------------
# Auxiliary functions for gwas analysis:
# --------------------------------------
domain = system("hostname -d", intern=TRUE)
if (length(domain) == 0) domain = ''
if (domain == 'broadinstitute.org'){
    bindir='~/data/EPIMAP_ANALYSIS/bin/'
} else {
    bindir='~/EPIMAP_ANALYSIS/bin/'
}
source(paste0(bindir, 'auxiliary_gwastree_plotting_functions.R'))
library(circlize)
library(dendextend)
library(ComplexHeatmap)
library(dplyr)
library(viridis)

# ----------------
# Get descendants:
# ----------------
get_dec <- function(subdend, declist, parent=0){
    node = attributes(unclass(subdend))$nodePar$pch
    # GET ID OF NODE:
    nset = declist$dec[[node]]
    # if (length(nset) == 0){ declist$dec[[node]] = NA }
    if (length(nset) == 0 || is.na(nset)){
        if (length(subdend) == 2){
            # Get subtrees and their node #s:
            dend1 = subdend[[1]]
            dend2 = subdend[[2]]
            d1 = attributes(unclass(dend1))$nodePar$pch
            d2 = attributes(unclass(dend2))$nodePar$pch
            # Update each node:
            declist = get_dec(dend1, declist, parent=node)
            declist = get_dec(dend2, declist, parent=node)
            # Merge each side's descendants:
            declist$dec[[node]] = c(declist$dec[[d1]], declist$dec[[d2]])
            declist$isleaf[node] = 0
        } else {
            # If leaf, get epigenome from matrix:
            slab = labels(subdend)
            declist$dec[[node]] = slab
            declist$isleaf[node] = 1
        }
        declist$parent[node] = parent
    } else { print("Already have consensus at this node") }
    return(declist)
}

# -------------------------
# Functions for enrichment:
# -------------------------
run.hyper <- function(y){phyper(q=y[1]-1, m=y[3], n=y[4] - y[3], k=y[2], lower.tail=FALSE)} 

run.fisher <- function(y){
    table <-  matrix(as.numeric(c(y[1], y[2], y[3], y[4])), ncol = 2, byrow = TRUE)
    if(any(is.na(table))) p <- "error" else p <- fisher.test(table, alternative="two.sided")$p.value
    return(p)
}

CUTP=12

get_pvalues = function(strait, dflist, cdlenlist, pdf, cutp=CUTP, 
                       ingroup='diff', outgroup='cons', against.parent=TRUE){
    # Find the GWAS uid with most snps:
    sgw = aggregate(pValue ~ uid + trait, gwdf[gwdf$trait == strait,], length)
    suid = sgw[order(sgw$pValue, decreasing=T), 'uid'][1]
    # Measure numbers of ingroup SNPs (in any group):
    isubdf = filter(dflist[[ingroup]], uid == suid)
    isnp = rep(0, NN)
    isnp[isubdf$node] = isubdf$nsnp
    ilen = cdlenlist[[ingroup]]
    # Make title:
    title = paste('Node', capitalize(ingroup),'vs.')
    # Measure numbers of outgroup SNPs (in any group + enh):
    if (outgroup == 'enh'){
        if (singlematch){
            # osnp = filter(nnumdf, uid == suid)$queryHits
            osnp = filter(nallnumdf, uid == suid)$queryHits
        } else {
            # osnp = filter(qnumdf, uid == suid)$queryHits
            osnp = filter(qallnumdf, uid == suid)$queryHits
        }
        olen = length(enhind)
        title = paste(title, 'All Enhancers')
    } else { 
        osubdf = filter(dflist[[outgroup]], uid == suid)  # SNPS in novel 
        osnp = rep(0, NN)
        osnp[osubdf$node] = osubdf$nsnp
        olen = cdlenlist[[outgroup]]
        if (against.parent){
            # Reorder snps and lengths by parent:
            osnp = osnp[declist$parent] 
            olen = olen[declist$parent]
            title = paste(title, 'Parent', capitalize(outgroup))
        } else {
            title = paste(title, 'Node', capitalize(outgroup))
        }
    }
    # Hyper-geometric testing:
    df = cbind(q=isnp, draw=ilen,
               m=osnp, N=olen)
    pout <- apply(df, 1, run.hyper)
    df = data.frame(df)
    df$pout = pout
    df$padj = p.adjust(df$pout)
    df$padj[df$padj == 0] = min(df$padj[df$padj != 0])
    df$rawlog10p = -log10(df$padj)
    # Cap and threshold the adjusted log10p:
    df$log10p = df$rawlog10p
    df$log10p[df$log10p < 3] = 0
    df$log10p[df$log10p > cutp] = cutp
    return(list(log10p=df$log10p, isnp=isnp, rawlp=df$rawlog10p, df=df, title=title))
}


# Test leaves against all enhancers:
get_pvalues_leaves = function(trait, dflist, cdlenlist, declist, cutp=CUTP){
    if (trait %in% gwdf$trait){
        sgw = aggregate(pValue ~ uid + trait, gwdf[gwdf$trait == trait,], length)
        suid = as.character(sgw[order(sgw$pValue, decreasing=T), 'uid'][1])
    } else if (trait %in% gwdf$uid){
        suid = trait
    } else {
        stop(paste0("Cannot find trait/uid:", trait))
    }
    # Measure numbers of consensus SNPs:
    isubdf = filter(dflist[['cons']], uid == suid)
    isnp = rep(0, NN)
    isnp[isubdf$node] = isubdf$nsnp
    ilen = cdlenlist[['cons']]
    # Measure numbers of enhancer snps:
    if (singlematch){
        osnp = filter(nallnumdf, uid == suid)$queryHits
    } else {
        osnp = filter(qallnumdf, uid == suid)$queryHits
    }
    olen = length(enhind)
    # Hyper-geometric testing:
    df = cbind(q=isnp, draw=ilen,
               m=osnp, N=olen)
    # Reduce to just leaves + test:
    ind = which(declist$isleaf == 1)
    df = df[ind,]
    pout <- apply(df, 1, run.hyper)
    df$node = ind  # Node number on tree
    df$pout = pout
    df$padj = p.adjust(df$pout)
    df$padj[df$padj == 0] = min(df$padj[df$padj != 0])
    df$rawlog10p = -log10(df$padj)
    # Cap and threshold the adjusted log10p:
    df$log10p = df$rawlog10p
    df$log10p[df$log10p < 3] = 0
    df$log10p[df$log10p > cutp] = cutp
    return(df)
}


# ------------------------------
# GWAS enrichment for a flat set
# ------------------------------
# INPUT: 
# trait
# dflist (all trait x flat interactions)
# lenlist
get_pvalues_flat = function(trait, dflist, lenlist, qdf, NF, cutp=CUTP, minp=3){
    # Traits:
    if (trait %in% qdf$uid){
        suid = trait
    } else {
        stop(paste0("Cannot find trait/uid:", trait))
    }
    # Measure numbers of consensus SNPs:
    isubdf = filter(dflist, uid == suid)
    isnp = rep(0, NF)
    isnp[isubdf$node] = isubdf$nsnp
    # Measure numbers of enhancer snps:
    osnp = filter(qallnumdf, uid == suid)$queryHits
    olen = length(enhind)
    # Hyper-geometric testing:
    df = cbind(q=isnp, draw=lenlist, m=osnp, N=olen)
    pout <- apply(df, 1, run.hyper)
    # Output dataframe:
    df = data.frame(node=1:NF, pout=pout)
    ll = pvalsdf.tolist(df, minp, cutp)
    ll$title = ""
    ll$isnp = isnp
    return(ll)
}


# Take enhancer sets + intersection dataframe
# Return number of intersecting snps per enhancer set:
eval.intersections <- function(enhsets, qdf){
    dflist = c()
    NL = length(enhsets)
    for (i in 1:NL){
        enhset = enhsets[[i]]
        # Map enhset to enhind (necessary for dflist)
        enhset = enhmap[enhset]
        enhtmp = merge(qdf, data.frame(subjectHits = enhset))
        if (nrow(enhtmp) > 0){
            # Count number of interactions (enh-snp) in set, add to dataframe:
            dfsnp = aggregate(queryHits ~ uid, enhtmp, length)
            names(dfsnp) = c('uid','nsnp')
            dfsnp$node = i
            dflist = rbind(dflist, dfsnp)
        }
    }
    return(dflist)
}


setup_dendrogram <- function(dend, ll, udf, declist, altline=3, cutp=CUTP, bcutoff=6, altlwd=c(1,1), palette=col3, ntop=5){
    log10p = ll$log10p
    # Color for strength of p-value:
    bins <- cut(log10p, seq(0, cutp, length.out=length(palette)), include.lowest=T) 
    cc = palette[bins]
    # Highlight the top hits:
    sord = order(ll$rawlp, ll$isnp, decreasing=T)
    tophits = head(sord, ntop)
    tophits = tophits[ll$log10p[tophits] > 0]  # Don't plot NS
    cexlist = .35 * (log10p > 0)
    cexlist[tophits] = .85
    pchlist = rep(19,NN)
    pchlist[tophits] = 21
    dend = set(dend, 'nodes_pch', pchlist)
    dend = set(dend, 'nodes_cex', cexlist)
    dend = set(dend, 'nodes_col', cc)
    # For counting:
    dend.dec = set(dend, 'nodes_pch', 1:NN)
    declist = list(dec=sapply(rep(NA, NN), list), isleaf=rep(NA, NN), parent=rep(NA, NN))
    declist = get_dec(dend.dec, declist=declist)
    # Get all descendant branches:
    keeplp = which((log10p >= bcutoff) * (1 - declist$isleaf) == 1) # Remove if "only" leaf
    hits = unique(as.character(unlist(sapply(keeplp, function(x){declist$dec[[x]]}))))
    # All hits and leaf hits: 
    leafhits = (log10p >= bcutoff)
    allhits = as.character(unlist(sapply(which(leafhits), function(x){declist$dec[[x]]})))
    if (length(allhits) > 0){
        hitdf = merge(data.frame(lab=allhits, nhit=1), data.frame(i=1:length(labels(dend)), lab=labels(dend)))
        hitdf = aggregate(nhit ~ lab + i, hitdf, sum)
    } else { hitdf = NULL }
    nhits = which(labels(dend) %in% hits)
    # NOTE: Could set altline to either 3 or 0 (blank - faster plotting)
    dend = set(dend, "by_labels_branches_lty", value=hits, TF_values = c(1, altline))
    dend = set(dend, "by_labels_branches_lwd", value=hits, TF_values = altlwd)
    # Edit colors:
    nl = which(declist$isleaf == 1)
    udf$col = 'grey65'
    udf$col[which(nl %in% which(leafhits == 1))] = 'red4'
    return(list(dend=dend, hits=leafhits, udf=udf, nhits=nhits, hitdf=hitdf))
}


# Is slow. can we speed up by not evaluating certain ones?
get_lr_pvalues = function(strait, dflist, cdlenlist, pdf, cutp=CUTP, 
                       ingroup='diff', outgroup='cons', against.parent=TRUE){
    # Trait and example place:
    sgw = aggregate(pValue ~ uid + trait, gwdf[gwdf$trait == strait,], length)
    suid = sgw[order(sgw$pValue, decreasing=T), 'uid'][1]
    NE = nrow(enhdf)
    pvals = rep(1, NN)
    for (i in 1:NN){
        print(paste(i, "is leaf?", declist$isleaf[i]))
        p = declist$parent[i]
        # Create vectors:
        allsnps = rep(0, NE)
        dann = rep(0, NE)
        pann = rep(0, NE)
        # Sets:
        if (singlematch){
            allsnpset = ndf[ndf$uid == suid, 'subjectHits']
        } else {
            allsnpset = qdf[qdf$uid == suid, 'subjectHits']
        }
        pset = enhmap[cdll[['cons']][[p]]]
        dset = enhmap[cdll[['diff']][[i]]]
        # Fill:
        dann[dset] = 1
        pann[pset] = 1
        allsnps[allsnpset] = 1
        # print(paste(sum(pann), sum(dann), sum(allsnps)))
        testdf = data.frame(snp=snps, parent=pann, node=ann, allsnp=allsnps, intercept=1)
        testdf$diff = testdf$node - testdf$parent
        # Log reg for each (~ 3.5 seconds)
        ptm <- proc.time()
        nmat = cbind(pann, 1)
        dmat = cbind(pann, dann, 1)
        nmodflr = fastLR(x=nmat, y=testdf$allsnp)
        dmodflr = fastLR(x=dmat, y=testdf$allsnp)
        lrdflr = dmodflr$loglikelihood - nmodflr$loglikelihood
        print(paste(sapply(c(nmodflr$loglikelihood, dmodflr$loglikelihood, lrdflr), digits=2, round)))
        pvals[i] = pchisq(lrdflr, df=1, lower.tail=FALSE)
        print(proc.time() - ptm)
    }
    # Make into dataframe:
    df = data.frame(node=1:NN, parent=declist$parent[1:NN])
    df$pout = pvals
    df$padj = p.adjust(df$pout)
    df$padj[df$padj == 0] = min(df$padj[df$padj != 0])
    df$rawlog10p = -log10(df$padj)
    # Cap and threshold the adjusted log10p:
    df$log10p = df$rawlog10p
    df$log10p[df$log10p < 3] = 0
    df$log10p[df$log10p > cutp] = cutp
    title = ''
    return(list(log10p=df$log10p, isnp=isnp, rawlp=df$rawlog10p, df=df, title=title))
}


# -------------------------------------------
# Fast implementation of logistic regression:
# -------------------------------------------
# Auxiliary:
sigm = function(z) {1 / (1 + exp(-z))}
rsigm = function(p) {log(p / (1 - p))}

calc_disjoint_lr = function(y, x1, x2=NULL, gamma=0.001, max.iter=1500, weights=NULL){
    NE = length(y)
    if (is.null(weights)){ weights = rep(1, NE) } 
    if (is.null(x2)){ x2 = rep(0, NE) }
    # Count statistics (could be weighted)
    N = length(y * weights)
    NY = sum(y * weights)
    NX1 = sum(x1 * weights)
    NTP1 = sum(x1 * y * weights)
    NFP1 = NX1 - NTP1 
    NX2 = sum(x2 * weights)
    NTP2 = sum(x2 * y * weights)
    NFP2 = NX2 - NTP2
    # For intercept, not the FP/TP:
    NI1 = NY - NTP1 - NTP2 # Y=1, X1=0, X2=0, I=1 (remaining snps)
    NI2 = N - NY - NFP1 - NFP2  # Y=0, X1=0, X2=0, I=1 (remaining loc - NY contains NTP1, NTP2)
    # Initialize coefficients:
    coeff = rep(0,3)
    pc = 1e-10 # with pseudo counts
    coeff[3] = rsigm((NY - NTP1 - NTP2)/(N - (NX1 + NX2)))
    coeff[1] = rsigm((NTP1 + pc)/(NX1)) - coeff[3]
    if (NX2 > 0){
        coeff[2] = rsigm((NTP2 + pc)/(NX2)) - coeff[3]
    } else {
        coeff[2] = 0
    }
    j = 1
    pctchange = 100
    # Update equations::
    while (pctchange > 0.0001 && j < max.iter){
        j = j + 1
        oldcoeff = coeff
        # Update x1 (1):
        inner1 = sigm(coeff[1] + coeff[3])
        inner2 = sigm(coeff[2] + coeff[3])
        inner3 = sigm(coeff[3])
        coeff[1] = oldcoeff[1] + 
            gamma * (NTP1 * (1 - inner1) + NFP1 * (0 - inner1))
        # Update x2 (2):
        if (NX2 > 0){
            coeff[2] = oldcoeff[2] +
                gamma * (NTP2 * (1 - inner2) + NFP2 * (0 - inner2))
        }
        # Update intercept (3):
        coeff[3] = oldcoeff[3] + 
            gamma * (NTP1 * (1 - inner1) + NFP1 * (0 - inner1) + 
                     NTP2 * (1 - inner2) + NFP2 * (0 - inner2) + 
                     NI1 * (1 - inner3) + NI2 * (0 - inner3))
        # Percent updated:
        pctchange = sum(abs(coeff - oldcoeff)) / sum(abs(oldcoeff))
    }
    # Calculate the log likelihood
    inner1 = sigm(coeff[1] + coeff[3])
    inner2 = sigm(coeff[2] + coeff[3])
    inner3 = sigm(coeff[3])
    llk = (NTP1 * log(inner1) + NFP1 * log(1 - inner1) + 
           NTP2 * log(inner2) + NFP2 * log(1 - inner2) + 
           NI1 * log(inner3) + NI2 * log(1 - inner3))
    # Return the parameters:
    return(list(loglikelihood=llk, coefficients=coeff))
}



# For where x2 is subset of x1
calc_subset_lr = function(y, x1, x2=NULL, gamma=0.001, max.iter=1500, weights=NULL){
    NE = length(y)
    if (is.null(weights)){ weights = rep(1, NE) } 
    if (is.null(x2)){ x2 = rep(0, NE) }
    # Count statistics (could be weighted)
    N = length(y * weights)
    NY = sum(y * weights)
    NX1 = sum(x1 * weights)
    NTP1 = sum(x1 * y * weights)  # Superset of NTP2
    NFP1 = NX1 - NTP1  # Superset of NFP2
    NX2 = sum(x2 * weights)
    NTP2 = sum(x2 * y * weights)
    NFP2 = NX2 - NTP2
    # For intercept, not the FP/TP:
    NI1 = NY - NTP1
    NI2 = N - (NY + NFP1)
    # Initialize coefficients:
    coeff = rep(0,3)
    pc = 0.000001  # with pseudo counts
    coeff[3] = rsigm(NY/N)
    coeff[1] = rsigm((NTP1 + pc)/(NTP1 + NFP1)) - coeff[3]
    if (NX2 > 0){
        coeff[2] = rsigm((NTP2 + pc)/(NTP2 + NFP2)) - coeff[3]
    } else {
        coeff[2] = 0
    }
    j = 1
    pctchange = 100
    # Update equations::
    while (pctchange > 0.0001 && j < max.iter){
        j = j + 1
        oldcoeff = coeff
        inner1 = sigm(coeff[1] + coeff[3])
        inner2 = sigm(coeff[1] + coeff[2] + coeff[3])  # Because x2 is subset of x1
        inner3 = sigm(coeff[3])
        # Update x1 (1):
        coeff[1] = oldcoeff[1] + gamma * ((NTP1 - NTP2) * (1 - inner1) + (NFP1 - NFP2) * (0 - inner1) +
                                          NTP2 * (1 - inner2) + NFP2 * (0 - inner2))
        # Update x2 (2):
        if (NX2 > 0){
            coeff[2] = oldcoeff[2] + gamma * (NTP2 * (1 - inner2) + NFP2 * (0 - inner2))
        }
        # Update intercept (3):
        coeff[3] = oldcoeff[3] + gamma * ((NTP1 - NTP2) * (1 - inner1) + (NFP1 - NFP2) * (0 - inner1) + 
                                       NTP2 * (1 - inner2) + NFP2 * (0 - inner2) + 
                                       NI1 * (1 - inner3) + NI2 * (0 - inner3))
        # Percent updated:
        pctchange = sum(abs(coeff - oldcoeff)) / sum(abs(oldcoeff))
    }
    
    # Calculate the log likelihood
    inner1 = sigm(coeff[1] + coeff[3])
    inner2 = sigm(coeff[1] + coeff[2] + coeff[3])
    inner3 = sigm(coeff[3])
    llk = ((NTP1 - NTP2) * log(inner1) + (NFP1 - NFP2) * log(1 - inner1) + 
           NTP2 * log(inner2) + NFP2 * log(1 - inner2) + 
           NI1 * log(inner3) + NI2 * log(1 - inner3))
    # Return the parameters:
    return(list(loglikelihood=llk, coefficients=coeff))
}


# For where x2 is subset of x1
calc_subset_lr = function(y, x1, x2=NULL, gamma=0.001, max.iter=1500, weights=NULL){
    NE = length(y)
    if (is.null(weights)){ weights = rep(1, NE) } 
    if (is.null(x2)){ x2 = rep(0, NE) }
    # Count statistics (could be weighted)
    N = length(y * weights)
    NY = sum(y * weights)
    NX1 = sum(x1 * weights)
    NTP1 = sum(x1 * y * weights)  # Superset of NTP2
    NFP1 = NX1 - NTP1  # Superset of NFP2
    NX2 = sum(x2 * weights)
    NTP2 = sum(x2 * y * weights)
    NFP2 = NX2 - NTP2
    # For intercept, not the FP/TP:
    NI1 = NY - NTP1
    NI2 = N - (NY + NFP1)
    # Initialize coefficients:
    coeff = rep(0,3)
    pc = 0.000001  # with pseudo counts
    coeff[3] = rsigm(NY/N)
    coeff[1] = rsigm((NTP1 + pc)/(NTP1 + NFP1)) - coeff[3]
    if (NX2 > 0){
        coeff[2] = rsigm((NTP2 + pc)/(NTP2 + NFP2)) - coeff[3]
    } else {
        coeff[2] = 0
    }
    j = 1
    pctchange = 100
    # Update equations::
    while (pctchange > 0.0001 && j < max.iter){
        j = j + 1
        oldcoeff = coeff
        inner1 = sigm(coeff[1] + coeff[3])
        inner2 = sigm(coeff[1] + coeff[2] + coeff[3])  # Because x2 is subset of x1
        inner3 = sigm(coeff[3])
        # Update x1 (1):
        coeff[1] = oldcoeff[1] + gamma * ((NTP1 - NTP2) * (1 - inner1) + (NFP1 - NFP2) * (0 - inner1) +
                                          NTP2 * (1 - inner2) + NFP2 * (0 - inner2))
        # Update x2 (2):
        if (NX2 > 0){
            coeff[2] = oldcoeff[2] + gamma * (NTP2 * (1 - inner2) + NFP2 * (0 - inner2))
        }
        # Update intercept (3):
        coeff[3] = oldcoeff[3] + gamma * ((NTP1 - NTP2) * (1 - inner1) + (NFP1 - NFP2) * (0 - inner1) + 
                                       NTP2 * (1 - inner2) + NFP2 * (0 - inner2) + 
                                       NI1 * (1 - inner3) + NI2 * (0 - inner3))
        # Percent updated:
        pctchange = sum(abs(coeff - oldcoeff)) / sum(abs(oldcoeff))
    }
    
    # Calculate the log likelihood
    inner1 = sigm(coeff[1] + coeff[3])
    inner2 = sigm(coeff[1] + coeff[2] + coeff[3])
    inner3 = sigm(coeff[3])
    llk = ((NTP1 - NTP2) * log(inner1) + (NFP1 - NFP2) * log(1 - inner1) + 
           NTP2 * log(inner2) + NFP2 * log(1 - inner2) + 
           NI1 * log(inner3) + NI2 * log(1 - inner3))
    # Return the parameters:
    return(list(loglikelihood=llk, coefficients=coeff))
}


# y = allsnps
# x = cbind(pann, dann)
# weights = NULL 

# General fast LR for observations and predictions in {0,1}
# Note: Here X is a matrix with ncol = nvar
calc_general_lr = function(y, x, gamma=0.001, max.iter=1500, weights=NULL){
    ptm <- proc.time()
    # Count statistics:
    N = length(y)
    NY = sum(y)
    NVAR = ncol(x)
    # if (is.null(weights)){ weights = rep(1, N) }
    if (NVAR == 1){ x = cbind(x, rep(0,N)) }
    # Get all contingency counts 
    # NOTE: May be slower than manually counting ~ .3s instead of .12:
    ctcounts = rep(0, 2^(NVAR + 1))
    for (vy in c(0, 1)){
        if (vy == 1){ iy = y } else { iy = 1 - y }
        for (v1 in c(0, 1)){
            if (v1 == 1){ i1 = x[,1] } else { i1 = 1 - x[,1] }
            for (v2 in c(0, 1)){
                if (v2 == 1){ i2 = x[,2] } else { i2 = 1 - x[,2] }
                # Store counts in binary (0:7 in 1:8):
                ind = (1 + vy * 1 + v1 * 2 + v2 * 4)
                if (is.null(weights)){
                    ctcounts[ind] = sum(iy * i1 * i2)
                } else {
                    ctcounts[ind] = sum(iy * i1 * i2 * weights)
                }
            }
        }
    }
    proc.time() - ptm
    # Initialize coefficients:
    coeff = rep(0,3)
    pc = 0.000001  # with pseudo counts
    coeff[3] = rsigm(NY/N)
    coeff[1] = rsigm((sum(ctcounts[1 + 1 + 2 + 4*0:1]) + pc) /
                     (sum(ctcounts[1 + 1 + 2 + 4*0:1]) + sum(ctcounts[1 + 0 + 2 + 4*0:1]))) - coeff[3]
    if (NVAR > 1){
        coeff[2] = rsigm((sum(ctcounts[1 + 1 + 2*0:1 + 4]) + pc) /
                         (sum(ctcounts[1 + 1 + 2*0:1 + 4]) + sum(ctcounts[1 + 0 + 2*0:1 + 4]))) - coeff[3]
    } else {
        coeff[2] = 0
    }
    j = 1
    pctchange = 100
    # Update equations::
    while (pctchange > 0.0001 && j < max.iter){
        j = j + 1
        oldcoeff = coeff
        # Update x1 (1):
        inner1 = sigm(coeff[1] + coeff[3])
        inner2 = sigm(coeff[2] + coeff[3])
        inner3 = sigm(coeff[3])
        coeff[1] = oldcoeff[1] + gamma * (NTP1 * (1 - inner1) + NFP1 * (0 - inner1))
        # Update x2 (2):
        if (NX2 > 0){
            coeff[2] = oldcoeff[2] + gamma * (NTP2 * (1 - inner2) + NFP2 * (0 - inner2))
        }
        # Update intercept (3):
        coeff[3] = oldcoeff[3] + gamma * (NTP1 * (1 - inner1) + NFP1 * (0 - inner1) + 
                                       NTP2 * (1 - inner2) + NFP2 * (0 - inner2) + 
                                       NI1 * (1 - inner3) + NI2 * (0 - inner3))
        # Percent updated:
        pctchange = sum(abs(coeff - oldcoeff)) / sum(abs(oldcoeff))
    }
    
    # Calculate the log likelihood
    inner1 = sigm(coeff[1] + coeff[3])
    inner2 = sigm(coeff[2] + coeff[3])
    inner3 = sigm(coeff[3])
    llk = (NTP1 * log(inner1) + NFP1 * log(1 - inner1) + 
           NTP2 * log(inner2) + NFP2 * log(1 - inner2) + 
           NI1 * log(inner3) + NI2 * log(1 - inner3))
    # Return the parameters:
    return(list(loglikelihood=llk, coefficients=coeff))
}


pvalsdf.tolist = function(df, minp, cutp){ 
    df$padj = p.adjust(df$pout)
    df$padj[df$padj == 0] = min(df$padj[df$padj != 0])
    df$rawlog10p = -log10(df$padj)
    # Cap and threshold the adjusted log10p:
    df$log10p = df$rawlog10p
    df$log10p[df$log10p < minp] = 0
    df$log10p[df$log10p > cutp] = cutp
    ll = list(log10p=df$log10p, rawlp=df$rawlog10p, df=df)
}


# ALTERNATE REGRESSION (slower):
# dmodel = calc_disjoint_lr(y=allsnps, x1=pann, x2=dann, gamma=gamma, weights=weights)
# pmat = cbind(parent=pann, intercept=1)
# dmat = cbind(parent=pann, diff=dann, intercept=1)
# dmat = cbind(diff=dann, pmat)
# NOTE: Can do cons-p vs. cons-p \\or// cons-p vs. cons-p + diff.
# Options tend to be similar.
# Log reg for each (~ 3.5 seconds)
# p2 = declist$parent[p]
# pmodflr = fastLR(x=pmat, y=y, start=pcoeff[[p2]])
# Better log reg for each (~ .01 seconds)

alt_lr = function(y, x1, x2=NULL, gamma=0.001, max.iter=1500, weights=NULL){
    require(RcppNumerical)
    # Make x matrix. (NOTE: y is allsnps vector)
    xmat = cbind(parent=x1, intercept=1)
    if (!is.null(x2)){ 
        xmat = cbind(parent=x1, diff=x2, intercept=1)
    }
    # Calculate model:
    model = fastLR(x=xmat, y=y)
    return(list(loglikelihood = model$loglikelihood,
                coefficients=model$coeff))
}


get_disjoint_lr = function(trait, type, against, qdf, minp=3, cutp=CUTP, 
                           weights=NULL, verbose=FALSE, return.coeff=FALSE, 
                           alt.lr=FALSE, max.iter=1500){
    # Trait and example place:
    if (trait %in% qdf$uid){
        suid = trait
    } else {
        if (trait %in% gwdf$trait){
            sgw = aggregate(pValue ~ uid + trait, gwdf[gwdf$trait == trait,], length)
            suid = sgw[order(sgw$pValue, decreasing=T), 'uid'][1]
        } else {
            stop(paste0("Cannot find trait/uid:", trait))
        }
    }
    NE = nrow(enhdf)
    pvals = rep(1, NN)
    pllk = rep(NA, NN)
    dllk = rep(NA, NN)
    diffllk = rep(NA, NN)
    isnp = rep(0, NN)
    if (return.coeff){ coeffmat = matrix(rep(0, NN * 3), nrow=NN, ncol=3) }
    pcoeff = sapply(rep(NA, NN), list)
    allsnpset = qdf[qdf$uid == suid, 'subjectHits']
    allsnps = rep(0, NE)
    allsnps[allsnpset] = 1
    title = paste('Node', capitalize(type),'vs.')
    if(against == 'enhancers'){
        title = paste(title, 'All Enhancers')
    } else {
        title = paste(title, capitalize(against), capitalize(type))
    }
    ptm <- proc.time()
    pb = txtProgressBar(min=0, max=NN, style = 3)
    for (i in 1:NN){
        setTxtProgressBar(pb, i)
        # Create vectors:
        nann = rep(0, NE)
        pann = rep(0, NE)
        # Sets:
        if (against == 'parent'){
            p = declist$parent[i]
            pset = enhmap[cdll[[type]][[p]]]
        } else if (against == 'sibling'){
            parent = declist$parent[i]
            d2 = which(declist$parent == parent)
            p = d2[d2 != i]
        } else if (against == 'enhancers'){
            pset = 1:NE
            p = 1 # Put the output here - same parent regression
            # All enhancers vs. all enhancers + node enhancers (separately)
            # - need to update calc_disjoint_lr to do NONDISJOINT.
        }
        nset = enhmap[cdll[[type]][[i]]]
        if (length(pset) == 0 || length(nset) == 0 || length(nset) == length(pset)){
            pvals[i] = 1
        } else {
            # Fill vectors:
            nann[nset] = 1
            pann[pset] = 1
            dann = nann - pann
            isnp[i] = sum(nann * allsnps)
            gamma = 0.01 # gamma = 0.00001
            if (type == 'union' || (!is.null(weights))){ gamma = gamma * 0.01 } # Tends to have high NA issues
            pll.condition = (is.na(pllk[p]) || pllk[p] > 0)
            while (is.na(pllk[p]) || pllk[p] > 0){ 
                if (alt.lr){
                    pmodel = alt_lr(y=allsnps, x1=pann)
                } else {
                    pmodel = calc_disjoint_lr(y=allsnps, x1=pann, gamma=gamma, weights=weights, max.iter=max.iter)
                }
                pcoeff[[i]] = pmodel$coefficients
                pllk[p] = pmodel$loglikelihood
                pll.condition = (is.na(pllk[p]) || pllk[p] > 0)
                if (pll.condition){
                    if (verbose) { print("NA parent LLK or bad regression, reduced gamma") } 
                    gamma = gamma * 0.1
                }
                if (gamma < 10^-10) break
            }
            gamma = 0.0001 
            if (type == 'union' || (!is.null(weights))){ gamma = gamma * 0.01 } # Tends to have high NA issues
            dll.condition = (is.na(dllk[i]) || is.infinite(dllk[i]) 
                             || dllk[i] > 0 || abs(diffllk[i]) > 5000)
            while (dll.condition){ 
                # Different purposes:
                if (type == 'union' || against == 'enhancers'){
                    if (alt.lr){
                        dmodel = alt_lr(y=allsnps, x1=pann, x2=nann)
                    } else {
                        dmodel = calc_subset_lr(y=allsnps, x1=pann, x2=nann, gamma=gamma, weights=weights)
                    }
                } else if (type == 'cons') {
                    if (alt.lr){
                        dmodel = alt_lr(y=allsnps, x1=pann, x2=dann)
                    } else {
                        dmodel = calc_disjoint_lr(y=allsnps, x1=pann, x2=dann, gamma=gamma, weights=weights, max.iter=max.iter)
                    }
                } else {
                    if (alt.lr){
                        dmodel = alt_lr(y=allsnps, x1=pann, x2=nann)
                    } else {
                        dmodel = calc_disjoint_lr(y=allsnps, x1=pann, x2=nann, gamma=gamma, weights=weights)
                    }
                }
                dllk[i] = dmodel$loglikelihood 
                diffllk[i] = 2 * (dllk[i] - pllk[p])  # NOTE: the LRT uses 2 x log(R1/R0)
                # Check condition:
                dll.condition = (is.na(dllk[i]) || is.infinite(dllk[i]) 
                                 || dllk[i] > 0 || abs(diffllk[i]) > 1000)
                if (dll.condition){
                    if (verbose) {print("NA diff LLK or bad regression, reduced gamma"); }
                    gamma = gamma * 0.1
                    if (verbose) {print(paste(sapply(c(i, pllk[p], dllk[i], diffllk[i]), digits=2, round)))}
                }
                if (gamma < 10^-10) break
            }
            pvals[i] = pchisq(diffllk[i], df=1, lower.tail=FALSE)
            if (diffllk[i] > 1000){
                print(paste0("Very large difference in log likelihood reported: Check: ", trait))
                print(paste(sapply(c(i, pllk[p], dllk[i], diffllk[i]), digits=2, round)))
            }
            if(return.coeff){ coeffmat[i,] = dmodel$coeff }
        }
    }
    close(pb)
    ptmend = proc.time() - ptm
    print(ptmend)
    # Make into dataframe:
    df = data.frame(node=1:NN, parent=declist$parent[1:NN], pout = pvals)
    ll = pvalsdf.tolist(df, minp, cutp)
    ll$title= title
    ll$isnp = isnp
    if (return.coeff){ ll$coeff = coeffmat }
    return(ll)
}


# -------------------------
# For adding genes to plot:
# -------------------------
packlist = function(x){
    # Make vertical list of two-three each
    x = sort(unique(x))
    str = ""
    # Pack, trying to keep to 4-5 lines max:
    ng = length(x)
    pack = 1
    if (ng > 8){ pack = 3 } else if (ng > 3){ pack = 2 }
    for (i in 1:(round(ng / pack - 1))){ 
        for (j in (pack*(i-1) + 1):(pack*i)){
            str = paste0(str, x[j])
            if (j == pack*i){
                str = paste0(str, '\n')
            } else {
                str = paste0(str, ', ')
            }
        }
    } 
    if (length(x) %% pack == 1){ str = paste0(str, x[length(x)]) }
    if (length(x) == 1){ str = x }
    return(str)
}

# Function for getting linked genes:
get.linked.genes = function(ll, suid, by.cons=TRUE, minp=3, allgenes=FALSE){
    # 1. Pull out the differential enhancers and their interactions: 
    # (associated with enrichment val)
    snpdf = qdf[qdf$uid == suid, ]
    names(snpdf) = c('snpid','enhid','uid')
    intset = sort(snpdf$enhid)
    redsnpdf = snpdf # For pruning
    ldf = ll$df
    ldf = ldf[order(ldf$pout),]
    if (!allgenes){
        ldf = ldf[ldf$rawlog10p > minp,]
    } 
    print(head(ldf,5))
    hitdf = c()
    for (j in 1:nrow(ldf)){
        if (nrow(redsnpdf) ==  0){ break } 
        i = ldf$node[j]
        lp = ldf$rawlog10p[j]
        p = declist$parent[i]
        pset = enhmap[cdll[[type]][[p]]]
        nset = enhmap[cdll[[type]][[i]]]
        dset = setdiff(nset, pset)
        iset = intersect(intset, dset)
        if (length(iset) > 0){
            idf = data.frame(enhid=iset, 
                             node=i,
                             log10p=lp,
                             uid=suid)
            # 2. Prune the SNPs - keep only top assoc per SNP (by pval):
            # Keeping only one SNP per node (but potentially multiple enhancers)
            idf = merge(idf, redsnpdf)
            # Remove the hit snps:
            snpset = unique(idf$snpid)
            redsnpdf = redsnpdf[!(redsnpdf$snpid %in% snpset),]
            intset = sort(redsnpdf$enhid)
            hitdf = rbind(hitdf, idf)
        }
    }
    nsnp = length(unique(snpdf$snpid))
    print(paste0(length(unique(hitdf$snpid))," of ",
                 nsnp, " snps assigned"))
    print(paste0(length(unique(hitdf$snpid[hitdf$log10p > minp]))," of ",
                 nsnp, " snps (interaction above minp = ", minp,")"))
    # 3. Get the locations of the enhancers:
    hitdf = cbind(hitdf, enhdf[hitdf$enhid, c('chr','start','end', 'cls')])
    hitdf$mid = (hitdf$start + hitdf$end)/2
    # Link to genes:
    if (!by.cons){
        # Option 1: Linked genes by linking agnostic to trees:
    } else {
        # Option 2: Linked genes by the RNA-seq consensus:
        # 3. Get consensus genes at each node:
        cgdf = c()
        for (i in unique(hitdf$node)){
            igenes = rnall$cons[[i]]
            itdf = tssdf[tssdf$gene %in% igenes, c('chr','tss', 'gene')]
            itdf$node = i
            cgdf = rbind(cgdf, itdf)
        }
        # 4. Match to nearest consensus gene at each node:
        cgdf = merge(cgdf, hitdf)
        cgdf$dist = abs(cgdf$tss - cgdf$mid)
        cgdf = merge(cgdf, aggregate(dist ~ node + enhid, cgdf, min))
    }
    # 5. Get the gene names per node + make plotting func:
    cgdf = merge(cgdf, gmdf)
    print(unique(cgdf$symbol))
    # Make legend:
    genetextdf = aggregate(symbol ~ node, cgdf, packlist)
    genetextdf = merge(genetextdf, nodedf)
    # Color (with the node:
    genetextdf = merge(genetextdf, unique(cgdf[,c('node','log10p')]))
    genetextdf$color = 'black'
    genetextdf$color[genetextdf$log10p < minp] = 'grey60'
    return(list(nodedf=genetextdf, cgdf=cgdf))
}


# Extract numbers (for sample-size):
munge.nos = function(x){
    locs = gregexpr("[0-9][0-9,]*",x)[[1]]
    ids = as.numeric(locs)
    lns = attributes(locs)$match.length
    # Extract:
    nos = sapply(1:length(ids), function(j){ 
                     num = substr(x, ids[j], ids[j]+ lns[j] - 1)
                     as.numeric(gsub(",","",num)) })
    return(nos)
}


# ------------------------------------------
# Method for cleaning up sample sizes:
# Associate each number with text - 
# look for case, control, individual, family
# Remove family and control.
# ------------------------------------------
prune.cases = function(x, only.cases=FALSE){
    locs = gregexpr("[0-9][0-9,]*",x)[[1]]
    ids = as.numeric(locs)
    lns = attributes(locs)$match.length
    # Extract:
    nos = sapply(1:length(ids), function(j){
                     num = substr(x, ids[j], ids[j]+ lns[j] - 1)
                     as.numeric(gsub(",","",num)) })
    # Extract text after each number:
    ids = c(ids, nchar(x) + 1)
    txt = sapply(1:(length(ids)-1), function(j){
                     substr(x, ids[j]+ lns[j], ids[j+1] - 1)
                     # as.numeric(gsub(",","",num)) 
                     })
    id.case = grep(" [Cc]ase", txt)
    id.indv = grep(" [Ii]ndividual", txt)
    id.control = grep(" [Cc]ontrol", txt)
    id.fam = grep(" [Ff]amil", txt)
    id.fam = id.fam[!(id.fam %in% id.case)]
    id.control = id.control[!(id.control %in% id.case)]
    id.rm = unique(c(id.fam, id.control))
    if (only.cases){
        if (length(id.case) > 0){
            nos = nos[id.case]
        } else { nos = 0 } 
    } else {
        if (length(id.rm) > 0){
            nos = nos[-id.rm]
        }
    }
    return(sum(nos))
}




# For controlling for FDR (against a permutation matrix):
fdr.adjust = function(ll, pmat, NPERM=300, NBELOW=3, cutp=3, maxp=12){
    # For 1% FDR
    pcut = apply(pmat[1:NPERM,], 2, function(x){sort(x)[NBELOW]})
    ldf = ll$df
    ldf$pcut = pcut
    # Adjust ll data:
    rlp = -log10(ldf$pout)
    nopass = !(ldf$pout < ldf$pcut)
    rlp[nopass] = 0
    lp = rlp
    # Cap non-raw values:
    lp[lp > maxp] = maxp
    lp[lp < cutp] = 0
    ll$rawlp = rlp
    ll$log10p = lp
    ll$df$rawlog10p = rlp
    ll$df$log10p = lp
    return(ll)
}

tsp.col = function(x, alpha=0.5){
    rr = col2rgb(x)
    rgb(t(rr)/255, alpha=alpha)
}

# Faster jaccard distance:
jacc.dist = function(x, y=NULL, verbose=TRUE){
    if (is.null(y)){ 
        if (verbose){ print("Using t(x) as second matrix") }
        y = t(x) 
    }
    rsx = apply(x, 1, sum)
    rsy = apply(y, 2, sum)
    NTX = length(rsx)
    NTY = length(rsy)
    rsx = matrix(rep(rsx, NTX), nrow=NTX,ncol=NTX, byrow=T)
    rsy = matrix(rep(rsy, NTY), nrow=NTY,ncol=NTY, byrow=T)
    rsm = t(rsx) + rsy
    nint = x %*% y
    nunion = rsm - nint
    dist = 1.0 - nint / nunion 
    return(dist)
}

# Cosine dist on two matrices:
cosine.dist = function(x, y=NULL, verbose=TRUE){
    if (is.null(y)){ 
        if (verbose){ print("Using t(x) as second matrix") }
        y = t(x) 
    }
    denom = sqrt(rowSums(x^2) %*% t(colSums(y^2)))
    dt = 1 - (x %*% y) / denom
    return(dt)
}

p.arrows2 = function(x1, y1, x2, y2, size = 1, width = (sqrt(5) - 1)/4/cin, fill = 2, border=NULL, ...){
    cin <- size * par("cin")[2]
    uin <- if (is.R()){ 1/xyinch()} else {par("uin")}
    segments(x1, y1, x2, y2, ...)
    x <- sqrt(seq(0, cin^2, length = floor(35 * cin) + 2))
    delta <- 0.005/2.54
    x.arr <- c(-x, -rev(x))
    wx2 <- width * x^2
    y.arr <- c(-wx2 - delta, rev(wx2) + delta)
    deg.arr <- c(atan2(y.arr, x.arr), NA)
    r.arr <- c(sqrt(x.arr^2 + y.arr^2), NA)
    theta <- atan2((y2 - y1) * uin[2], (x2 - x1) * uin[1])
    lx <- length(x1)
    Rep <- rep.int(length(deg.arr), lx)
    x2 <- rep.int(x2, Rep)
    y2 <- rep.int(y2, Rep)
    theta <- rep.int(theta, Rep) + rep.int(deg.arr, lx)
    r.arr <- rep.int(r.arr, lx)
    # Arrowheads:
    polygon(x2 + r.arr * cos(theta)/uin[1], y2 + r.arr * sin(theta)/uin[2],
            col = fill, border=border)
}
