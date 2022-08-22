#!/usr/bin/R
# ---------------------------------------------------
# Preliminary analysis for tree GWAS AND linked genes 
# Counts on a fixed matrix from epigenomic similarity
# RNA-seq from ENCODE data
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

gtdir = "gwas_tree_analysis/"
regdir = paste0(gtdir, "regressions/")
perdir = paste0(gtdir, "permuted_catalogs/")
epref = paste0(usetree, '_e', tol, '_')
regpref = paste0(regdir, epref)
perpref = paste0(perdir, epref)
cmd = paste('mkdir -p', gtdir, regdir, perdir)
system(cmd)

# ------------------
# Load RNA-seq data:
# ------------------
# Read in the RNA-seq protein-coding genes:
rnadir = 'RNA-seq/files/RNA-seq/tsv/'
lmatpcfile = paste0(rnadir, 'merged_log2fpkm.pc.mtx.gz')
availfile = paste0(rnadir, 'avail_rna.txt')
if (!file.exists(availfile)){
    pcdf = read.delim(gzfile(lmatpcfile), header=T, sep="\t")
    pcid = unique(pcdf$id)
    write.table(pcid, availfile, quote=F, row.names=F, sep="\t")
} else { 
    pcid = scan(availfile, 'c', sep="\n")
}
rna.avail = 1 * (labels(dend) %in% pcid)

# Load in TSS annotation:
tssfile = 'Annotation/gencode.gene.tss.v30lift37.basic.annotation.gtf.gz'
tssdf = read.delim(gzfile(tssfile), header=F, sep="\t", stringsAsFactors=F)
names(tssdf) = c('chr','tss','t2','gene')

gmfile = 'Annotation/gencode.gene.name.mapping.tsv'
gmdf = read.delim(gmfile, header=F, sep="\t", stringsAsFactors=F)
names(gmdf) = c('gene', 'symbol')

# --------------------------------------------------------
# Make RNA-seq consensus at a variety of log2fpkm cutoffs:
# --------------------------------------------------------
rnallfile = paste0(rnadir, 'rna_consensus_object_', usetree, '_070819.Rdata')
if (!file.exists(rnallfile)){
    # Recursive function to get consensus and merge:
    subpcdf = pcdf[pcdf$log2fpkm > 3,]
    subpcgenes = unique(pcdf$gene)
    pcid = unique(pcdf$id)
    get_consensus_rna <- function(subdend, rnall){
        node = attributes(unclass(subdend))$nodePar$pch
        # Get id of node:
        print(node)
        nset = rnall$cons[[node]]
        print(head(nset))
        if (length(nset) == 0 || is.na(nset)){
            if (length(subdend) == 2){
                # If internal node, get consensus by merging:
                # Get subtrees and their node #s:
                dend1 = subdend[[1]]
                dend2 = subdend[[2]]
                d1 = attributes(unclass(dend1))$nodePar$pch
                d2 = attributes(unclass(dend2))$nodePar$pch
                # Update each node:
                rnall = get_consensus_rna(dend1, rnall)
                rnall = get_consensus_rna(dend2, rnall)
                # Merge nodes for consensus (intersect) or union:
                rnall$cons[[node]] = intersect(rnall$cons[[d1]], rnall$cons[[d2]])
                rnall$union[[node]] = union(rnall$union[[d1]], rnall$union[[d2]])
                # Update each of the descendants - diff is in node, not parent
                rnall$diff[[d1]] = setdiff(rnall$cons[[d1]], rnall$cons[[node]])
                rnall$diff[[d2]] = setdiff(rnall$cons[[d2]], rnall$cons[[node]])
                # Update: novel is in parent, not node
                rnall$novel[[d1]] = setdiff(rnall$union[[node]], rnall$union[[d1]])
                rnall$novel[[d2]] = setdiff(rnall$union[[node]], rnall$union[[d1]])
            } else {
                print('leaf')
                # If leaf, get epigenome from matrix:
                id = labels(subdend)
                if (id %in% pcid){
                    rnall$cons[[node]] = subpcdf$gene[subpcdf$id == id]
                } else {
                    rnall$cons[[node]] = subpcgenes
                }
                rnall$union[[node]] = rnall$cons[[node]]
            }
        } else { print("Already have consensus at this node") }
        print(length(rnall$cons[[node]]))
        return(rnall)
    }
    # Building from the bottom, fill in all:
    rnall = list(cons=sapply(rep(NA, NN), list),
                 diff=sapply(rep(NA, NN), list),
                 union=sapply(rep(NA, NN), list),
                 novel=sapply(rep(NA, NN), list))
    # Store node id as pch (use to ID where we are)
    set(dend, 'nodes_pch', 1:NN) -> dend
    # Run consensus function to update clist:
    rnall = get_consensus_rna(dend, rnall)
    save(rnall, file=rnallfile)
} else {
    print("[STATUS] Loading rnall from file")
    load(rnallfile)
}

# NOTE: Will need to try to fill a lot of gaps on RNA-seq.
# TODO: Plotting track showing if had RNA-seq or not


# ---------------------------------------
# Look at specific GWAS to pilot linking:
# ---------------------------------------
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
# suid = "29187748 - Lateral ventricle volume in trauma-exposed individuals"
suid = "26343387 - Myocardial infarction"

suid = "30617256 - Alzheimer's disease or family history of Alzheimer's disease"

# Logistic regression arguments:
# type = 'union'
# type = 'novel'
# type = 'diff'
type = 'cons'
against = 'parent'
weighted = FALSE
apref = paste0(type, '_', against)
if (weighted){ 
    weights = sqrt(1 / matmarg[,2])
    apref = paste0(apref, '_weighted')
} else {
    weights = NULL
}

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

nsnpdf = aggregate(pValue ~ uid, gwdf, length)
NCAT = 300

# Run or load regression:
pdf(paste0(treeimgpref,apref, '_examples_small.pdf'), width=3.5,height=4, onefile=T)
MINP=3
for (suid in suidlist) {
    plotind = which(uids == suid)
    ipref = paste0(apref, sprintf("_%05d", plotind))
    print(paste(suid, '----- output with prefix', ipref))

    # Perform or load regression:
    suffix = '_adj300_1.Rda'
    cregfile = paste0(regpref, ipref, '_lreg', suffix)
    # regfile = paste0(regpref, ipref, '_lreg.Rda')
    # TODO: Remove me: just for speed
    if (file.exists(cregfile)){
        # Perform or load regression:
        regfile = paste0(regpref, ipref, '_lreg.Rda')
        if (!file.exists(regfile)){
            ll = try(get_disjoint_lr(trait=suid, qdf=qdf, type=type, against=against, weights=weights))
            if (class(ll) == 'try-error'){
                ll = 'error'
            }
            save(ll, file=regfile)
        } else {
            load(regfile)
        }

        # -----------------------
        # Correct the regression:
        # -----------------------
        if (ll != 'error'){
            # Load in the null p-values matrix:
            nsnps = nsnpdf$pValue[nsnpdf$uid == suid]
            print("Adjusting to 0.3%")
            npref = paste0(apref, sprintf("_permuted_%06d", nsnps))
            pmatfile = paste0(perpref, npref, '_', NCAT, '_lreg_pmat.Rda')
            load(pmatfile)
            # Try p < 0.3% FDR
            ll = fdr.adjust(ll, pmat, NPERM=NCAT, NBELOW=1)
        }

        # If plot at diff minp:
        cutp = 12
        if (MINP != 3){
            df = ll$df
            df$log10p = df$rawlog10p
            df$log10p[df$log10p < MINP] = 0
            df$log10p[df$log10p > cutp] = cutp
            ll$log10p=df$log10p
            ll$rawlp=df$rawlog10p
            ll$df = df
        }

        # Add the sample size:
        init = unique(gwssdf$initSample[gwssdf$uid == suid])
        repl = unique(gwssdf$replSample[gwssdf$uid == suid])
        strait = unique(gwdf$trait[gwdf$uid == suid])
        # Make sample size onto two lines:
        init = strsplit(init, ", ")[[1]]
        md = round((length(init) + 1) /2)
        init = paste0(paste0(init[1:md], collapse=", "), "\n",
                      paste0(init[md:length(init)], collapse=", "))

        # Link genes (optional), setup dendrogram:
        # gg = get.linked.genes(ll, suid, allgenes=FALSE, minp=MINP)
        dd  = setup_dendrogram(dend3, ll, udf, declist=declist, bcutoff=MINP, altline=0)

        # ---------------------
        # Plot small dendrogram
        # ---------------------
        NTOP=5
        ldf = ll$df
        ldf = ldf[order(ldf$pout),]
        ldf = ldf[1:NTOP,] # TOP N only
        ldf = ldf[ldf$rawlog10p > MINP,]
        ntdf = merge(ldf, nodetissue)
        ntdf = ntdf[,c('node','GROUP','COLOR')]
        names(ntdf) = c('node','symbol','color')
        ntdf = merge(ntdf, nodedf)
        ntdf$cex=.5

        par(mar=c(1,0,1,0))
        circleplot(dend=dd$dend, lab=lab, udf=dd$udf, fractional=FALSE, hit.track=dd$hitdf, 
                   plot.labels=FALSE, add.metadata=TRUE, reduced.metadata=TRUE, nodedf=ntdf, scale=2)
        circos.clear()
        mtext(strait, side=3, line=0, cex=.7)
        mtext(init, side=1, line=0, cex=.25)
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
# circleplot(dend=dd$dend, lab=lab, udf=dd$udf, fractional=FALSE, hit.track=dd$hitdf, leaf.pval=leafmat[strait,], nodedf=gg$nodedf, rna.avail=rna.avail)
circleplot(dend=dd$dend, lab=lab, udf=dd$udf, fractional=FALSE, hit.track=dd$hitdf, leaf.pval=leafmat[strait,])
upViewport()
circos.clear()
# draw(pd.legend.ext, x = circle_size, just = "left")
title(paste(suid, '-', ll$title))
dev.off()



