#!/usr/bin/R
# ------------------------------------------------------
# Plot the SNPs and the GO terms (+ small figure)
# For ALL 538 traits. (+ recalculate)
# + for ALL 836 traits at 1% as well (show top 20 nodes)
# qsub -cwd -t 1-55 -P compbio_lab -l h_vmem=25G -l h_rt=3:00:00 -hold_jid flat_gtchunked -N flat_calcgtchunked -j y -b y -V -r y -o $DBDIR/out/ChromImpute "source activate mv_env; R --slave -f $BINDIR/calculate_all_flat_enrichments.R"
# ------------------------------------------------------
domain = system("hostname -d", intern=TRUE)
if (length(domain) == 0) domain = ''
if (domain == 'broadinstitute.org'){
    bindir='~/data/EPIMAP_ANALYSIS/bin/'
} else {
    bindir='~/EPIMAP_ANALYSIS/bin/'
}
source(paste0(bindir, 'general_EPIMAP_ANALYSIS.R'))
# library(qvalue)
# library(flexclust)
# library(gridExtra)
# library(kableExtra)  # For latex table.
# library(grid)
library(GenomicRanges)

# # ----------------------
# # Arguments for modules:
# # ----------------------
# filepref = 'cls_merge2_wH3K27ac100_300'
# tagline = 'ChromHMM Enhancers'
# imgdir = paste0(img, "clusters/") 
# # Load in and process data (saves to matrices):
# commandArgs <- function(trailingOnly=TRUE){
#     c(filepref, tagline, imgdir) }
# source(paste0(bindir, 'load_modules_data.R'))

# ------------------------
# Arguments for GWAS data:
# ------------------------
usetree = 'enhancers'
tol = 2500
singlematch = FALSE
plotting.only = FALSE  # No, need to run the regressions
use.adj = TRUE
use.strict = TRUE
MINP=3

# Load in + process all of the relevant matrices/datasets:
commandArgs <- function(trailingOnly=TRUE){
    c(usetree, tol, singlematch, plotting.only, use.adj, use.strict) }
source(paste0(bindir, 'load_statistics_gwastree_enrichments.R'))

# ------------
# Directories:
# ------------
gtdir = "gwas_tree_analysis/"
regdir = paste0(gtdir, "regressions/")
perdir = paste0(gtdir, "permuted_catalogs/")
caddir = paste0(img, 'CAD_example/')
epref = paste0(usetree, '_e', tol, '_')
cadpref = paste0(caddir, 'ALLTRAITS_')
regpref = paste0(regdir, epref)
perpref = paste0(perdir, epref)
cmd = paste('mkdir -p', gtdir, regdir, perdir, caddir)
system(cmd)

# ------------------
# Load RNA-seq data:
# ------------------
# Read in the RNA-seq protein-coding genes:
rnadir = 'RNA-seq/files/RNA-seq/tsv/'
lmatpcfile = paste0(rnadir, 'merged_log2fpkm.pc.mtx.gz')
availfile = paste0(rnadir, 'avail_rna.txt')
if (!file.exists(availfile)){
    pcdf = read.delim(gzfile(lmatpcfile), header=T, sep="\t", stringsAsFactors=F)
    pcid = unique(pcdf$id)
    write.table(pcid, availfile, quote=F, row.names=F, sep="\t")
} else { 
    pcid = scan(availfile, 'c', sep="\n")
}
rna.avail = 1 * (labels(dend) %in% pcid)

# (very rough) average per tissue:
avgfile = paste0(rnadir, 'avg_rna.Rda')
if (!file.exists(avgfile)){
    pcdf = merge(pcdf, meta[,c('id','GROUP')])
    red.pcdf = aggregate(log2fpkm ~ gene + GROUP, pcdf, function(x){x = 2^x; log2(mean(x))})
    rm(pcdf)
    save(red.pcdf, file=avgfile)
} else {
    load(avgfile)
}

# Load in TSS annotation:
tssfile = 'Annotation/gencode.gene.tss.v30lift37.basic.annotation.gtf.gz'
tssdf = read.delim(gzfile(tssfile), header=F, sep="\t", stringsAsFactors=F)
tssgr = GRanges(tssdf$chr, IRanges(tssdf$tss, tssdf$tss))
names(tssdf) = c('chr','tss','t2','gene')

gmfile = 'Annotation/gencode.gene.name.mapping.tsv'
gmdf = read.delim(gmfile, header=F, sep="\t", stringsAsFactors=F)
names(gmdf) = c('gene', 'symbol')

gwrsdf = read.delim(gwrsidfile, header=T, sep="\t", stringsAsFactors=F)
gwrsdf = aggregate(name ~ chrom + chromStart + uid, gwrsdf, function(x){paste(sort(unique(x)), collapse=", ")})
gwdf2 = merge(gwdf, gwrsdf, all.x=TRUE)

# -------------------------------------
# Load all gwas enrichment information:
# -------------------------------------
load(paste0(eprefix, 'kept_allgwas_ordered', suffix))
print(length(Znam))  # Traits to look at/plot

snpfile = paste0(regdir, eprefix, apref, '_logreg_all_wsnp', suffix)
load(snpfile)
nodetissue = nodetissue[order(nodetissue$node),]
MINP=3

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
leafrep[lind] = nodetissue$GROUP[lind]


# Run analysis for all (no plotting yet):
# suid = Znam[513]
for (suid in Znam[89:538]){
    print(suid)
    print(which(Znam == suid))
    suid.lp = rmat[suid,]
    sub.qdf = qdf[qdf$uid == suid, ]
    # Pick up top 20 nodes + their metadata:
    top.nodes = which(suid.lp > 3)
    NPLOT = min(20, length(top.nodes))  # TOP N
    nodegroup = nodetissue$GROUP[top.nodes]
    nodecol = nodetissue$COLOR[top.nodes]
    rawnodelp = suid.lp[top.nodes]
    nodelp = round(rawnodelp,1)
    # nodenames = paste0(nodegroup, ' - ', nodelp)
    nodenames = paste0(top.nodes, '. ', nodegroup, ' - ', nodelp)
    nord = tail(order(suid.lp[top.nodes]), NPLOT)
    plotind= which(uids == suid)
    ipref = paste0(apref, sprintf("_%05d", plotind))
    suid.leafrep = leafrep[top.nodes]

    # ----------------------
    # Make the SNP matrices:
    # ----------------------
    if (length(top.nodes) > 1){
        out = sapply(top.nodes, function(x){cdll$diff[[x]]})
        omat = sapply(out, function(x){ x = enhmap[x]; return(1 * (sub.qdf$subjectHits %in% x)) }) 
    } else {
        out = enhmap[cdll$diff[[top.nodes]]]
        omat = as.matrix(1 * (sub.qdf$subjectHits %in% out))
    }
    colnames(omat) = nodenames
    tform = make.tform(sub.qdf$queryHits)
    smat = t(tform) %*% omat
    smat[smat > 0] = 1
    snps = unique(sub.qdf$queryHits)
    rownames(smat) = paste0('s', snps)

    set.seed(1)
    smat = smat[apply(smat, 1, sum) > 0,, drop=F]
    if (ncol(smat) > 1 && nrow(smat) > 2){
        dnj.cons = dist(smat, 'jaccard')
        hnj.cons <- hclust(dnj.cons, method='ward.D')
        cocl = order.optimal(as.dist(dnj.cons), hnj.cons$merge)$order
        snp.cons.rn <- names(cocl)[cocl]
        snp.cons.rn = snp.cons.rn[snp.cons.rn %in% rownames(smat)]
        NCUT = ifelse(nrow(smat) < 20, min(nrow(smat)-2,5), 10)
        snp.cto = cutree(hnj.cons, NCUT)[cocl]
    } else {
        snp.cons.rn = rownames(smat)
        snp.cto = rep(1, nrow(smat))
        names(snp.cto) = snp.cons.rn
    }

    # ----------------------
    # GO terms for clusters:
    # ----------------------
    snpdf = data.frame(queryHits=as.numeric(sub("s","",names(snp.cto))), 
                       snp=names(snp.cto), cls=snp.cto)
    snpdf$chr = paste0('chr', gwdf$chrom[snpdf$queryHits])
    snpdf$start = gwdf$chromStart[snpdf$queryHits]
    snpdf$rsid = gwdf2$name[snpdf$queryHits]
    snpgr = GRanges(snpdf$chr, IRanges(snpdf$start, snpdf$start))
    snpdf$gene = tssdf$gene[nearest(snpgr, tssgr)]
    stssdf = merge(snpdf, gmdf)
    # Aggregate GO test sets per cluster:
    tree.setsdf = aggregate(symbol ~ cls, stssdf,
                            function(x){paste(unique(x), collapse=", ")})
    mapsnpdf = aggregate(cbind(symbol, rsid) ~ snp, stssdf, function(x)x[1])
    mapsnp = paste0(mapsnpdf$symbol, ' (', mapsnpdf$rsid, ')')
    names(mapsnp) = mapsnpdf$snp

    plot_snp_full = function(mat, numsnp=NULL, leafrep=leafrep){
        sp=0.1
        bsp=3
        tsp=10
        # Side 1
        xat = seq(0,1,length.out=NPLOT)
        par(mar=c(sp, bsp, tsp, sp))
        image(t(mat), col=c('white','grey25'), axes=F, useRaster=T)
        text(y=par()$usr[4] + .005 * diff(par()$usr[3:4]), x=xat, 
             paste0(1:length(leafrep), '. ', leafrep), col=rev(nodecol[nord]), 
             srt=90, adj=0, xpd=TRUE, cex=.5)
        box(lwd=.5)
        msnp = mapsnp[rownames(mat)]
        if (!is.null(numsnp)){ msnp = paste0(rev(numsnp) , '. ', msnp) }
        yat = seq(0, 1, length.out=length(msnp))
        if (length(msnp) > 200) {txtcex=.1} else {txtcex=.3}
        text(x=par()$usr[1]-.01 * diff(par()$usr[1:2]), y=yat,
             msnp, col='black', srt=0, adj=1, xpd=TRUE, cex=txtcex)
    }

    # pdf(paste0(cadpref, ipref, 'all_snps_top20.pdf'), width=2, height=6)
    # par(xaxs='i')
    # par(yaxs='i')
    # sp=.15
    # smat.reord = smat[snp.cons.rn, rev(nodenames[nord])]
    # smat.reord = smat.reord[apply(smat.reord,1, sum) > 0,]
    # plot_snp_full(smat.reord, numsnp=1:nrow(smat.reord), leafrep=rev(leafrep[top.nodes[nord]]))
    # dev.off()

    # --------------------------------
    # Aggregate GO test sets per node:
    # --------------------------------
    print("[STATUS] Getting nearest expressed gene")
    smatdf = data.frame(t(smat[,nodenames[nord], drop=F]))
    smatdf$node = nodenames[nord]
    smatdf = gather(smatdf, snp, value, -node)
    smatdf = smatdf[smatdf$value > 0,]
    smatdf$pval = as.numeric(sub(".* - ", "",smatdf$node))
    smatdf$GROUP = sub(" - .*", "",smatdf$node)
    smarg = apply(smat, 1, sum)
    smargdf = data.frame(snp=names(smarg), marg=smarg)
    smatdf = merge(smatdf, smargdf)
    smatdf$queryHits = as.numeric(sub("s","",smatdf$snp))
    smatdf$chr = paste0('chr',gwdf$chrom[smatdf$queryHits])
    smatdf$start = gwdf$chromStart[smatdf$queryHits]

    # Add gene expr, tss (this step usually rate-limiting):
    stssdf = merge(tssdf, smatdf) # ALL.
    stssdf$dist = abs(stssdf$tss - stssdf$start)
    # hist(red.pcdf$log2fpkm, 50)
    red.pcdf = red.pcdf[red.pcdf$log2fpkm > 2,]
    stssdf = merge(stssdf, red.pcdf)
    # Add nearest (expressed gene)
    if (nrow(stssdf) > 1){
        stssdf = merge(stssdf, aggregate(dist ~ snp, stssdf, min))
    }
    stssdf = merge(stssdf, gmdf)
    # If any snp missing, add the nearest gene:
    smatdf = merge(smatdf, mapsnpdf)
    names(smatdf)[ncol(smatdf)] = 'nearest.symbol'
    stssdf = merge(smatdf, stssdf, all.x=TRUE)
    sind = is.na(stssdf$symbol)
    stssdf$symbol[sind] = stssdf$nearest.symbol[sind]
    # Node centered sets:
    tree.nodesdf = aggregate(symbol ~ node + pval, stssdf, function(x){paste(unique(x), collapse=", ")})
    tree.nodesdf = tree.nodesdf[order(tree.nodesdf$pval, decreasing=T),]

    print("[STATUS] Running GO on nodes")
    tn.nodefile = paste0(regpref, ipref, '_go_nodes', suffix)
    tn.resdf = c()
    fulltn.resdf = c()
    if (!file.exists(tn.nodefile)){
        for (i in 1:nrow(tree.nodesdf)){
            x = tree.nodesdf$symbol[i]
            x = strsplit(x, ", ")[[1]]
            df = try(go.enr(x, minsize=10, maxsize=150))
            if (class(df) == 'try.error'){ df = NULL }
            if (!is.null(df)){
                node = tree.nodesdf$node[i]
                df$node = node
                df$GO = rownames(df)
                fulltn.resdf = rbind(fulltn.resdf, df)
                tn.resdf = rbind(tn.resdf, df[1:min(3,nrow(df)),])
            }
        }
        if (!is.null(tn.resdf)){
            tn.resdf$geneID = gsub("/",", ", tn.resdf$geneID)
            tn.resdf$geneID = sapply(tn.resdf$geneID, width=50, split.text)
            tn.resdf$Description = sapply(tn.resdf$Description, width=90, split.text)
            tn.resdf$pval = as.numeric(sub("^.* - ", "", tn.resdf$node))
            tn.resdf = tn.resdf[order(tn.resdf$pval, decreasing=T), ] 
        }
        save(fulltn.resdf, tn.resdf, file=tn.nodefile)
    } else {
        load(tn.nodefile)
    }

    # Prune terms by qvalue first:
    if (!is.null(fulltn.resdf)){
        sub.godf = fulltn.resdf[, c('node','Description', 'qvalue', 'GO', 'geneID')]
        sub.godf = sub.godf[!is.na(sub.godf$qvalue),]
    }
    if (!(is.null(sub.godf)) && nrow(sub.godf) > 0){
        if (nrow(sub.godf) > 200 || min(sub.godf$qvalue, na.rm=T) < 1e-4){
            sub.godf = sub.godf[sub.godf$qvalue < 1e-2,]
        } else {
            sub.godf = sub.godf[sub.godf$qvalue < 1e-1,]
        }

        sub.godf = sub.godf[nchar(sub.godf$Description) < 100,]
        ndesc = aggregate(node ~ GO + Description, sub.godf, length)
        ndesc = ndesc[order(ndesc$node, decreasing=T), ] 
        go.cutoff = 10
        go.desc = ndesc$Description[ndesc$node <= go.cutoff]
        go.id = ndesc$GO[ndesc$node <= go.cutoff]

        # Keep 1-3 terms per term, make dataframe with all incidences:
        sub.godf = sub.godf[order(sub.godf$qvalue),]
        sub.godf.old = sub.godf
        sub.godf = sub.godf[sub.godf$Description %in% go.desc,]
        # go.keep = go.desc
        go.keep = c()
        NKEEP = 6
        for(ni in nodenames[nord]){ 
            go.keep = c(go.keep, head(sub.godf$Description[sub.godf$node == ni], NKEEP)) }
        go.keep = unique(go.keep)

        process_godf = function(ssdf){
            ssdf$log10q = -log10(ssdf$qvalue)
            sswide = spread(ssdf[,c('Description','node','log10q')], node, log10q, fill=0)
            mat = as.matrix(sswide[,-1, drop=F])
            rownames(mat) = sswide$Description
            rn = as.character(rev(nodenames[nord]))
            rn = rn[rn %in% colnames(mat)]
            mat = mat[, rn, drop=F]
            if (ncol(mat) > 1){
                rd = diag.mat2(t(mat) > 0)
                go.ord = rd[[2]]
            } else {
                go.ord = rownames(mat)
            }
            # Make plot:
            ssdf$node = factor(ssdf$node, levels=nodenames[nord])
            ssdf$Description = factor(ssdf$Description, levels=rev(go.ord))
            ssdf$cex = 1.5
            ssdf$cex[ssdf$qvalue < 1e-4] = 1
            ssdf$cex[ssdf$qvalue < 1e-3] = .75
            lindf = rbind(data.frame(x0=1:length(go.ord), x1=1:length(go.ord), y0=1, y1=NPLOT), 
                          data.frame(x0=1, x1=length(go.ord), y0=1:NPLOT, y1=1:NPLOT))
            ylindf = data.frame(x0=1:length(go.ord), x1=1:length(go.ord), y0=1, y1=NPLOT)
            return(list(ssdf, lindf, ylindf, go.ord))
        }

        # Get order for GO terms:
        sub.godf = sub.godf[sub.godf$Description %in% go.keep,]
        ll = process_godf(sub.godf)
        ssdf = ll[[1]]
        lindf = ll[[2]]
        ylindf = ll[[3]]
        go.ord = ll[[4]]

        # ---------------------------------------
        # Plot the GO terms as a separate figure:
        # ---------------------------------------
        plot_go_full = function(rnrow, leafrep=leafrep){
            sp=0.1
            bsp=12
            tsp=8
            go.genedf = aggregate(geneID ~ Description, ssdf, function(x){x[which.max(nchar(x))]})
            rownames(go.genedf) = go.genedf$Description
            go.genedf = go.genedf[go.ord,]
            ssdf$row = length(go.ord) - as.numeric(ssdf$Description) + 1 
            sub.ssdf = ssdf[ssdf$row %in% rnrow,]
            go.ord.sub = go.ord[rnrow]
            lindf = rbind(data.frame(x0=rnrow, x1=rnrow, y0=1, y1=NPLOT), 
                          data.frame(x0=min(rnrow), x1=max(rnrow), y0=1:NPLOT, y1=1:NPLOT))
            par(mar=c(sp, bsp, tsp, sp))
            ylim=c(min(rnrow) - .5, max(rnrow) + .5)
            plot(1, type='n', axes=F, ylab='', xlab='', 
                 xlim=c(.5,NPLOT+.5), ylim=ylim)
            segments(y0=lindf$x0, y1=lindf$x1, 
                     x0=lindf$y0, x1=lindf$y1, 
                     lwd=.1, lty='solid')
            text(y=rnrow, 
                 x=par()$usr[1] - 0.01 * diff(par()$usr[1:2]),
                 paste0(go.ord.sub, '\n(', go.genedf[go.ord.sub,'geneID'], ')'),
                 srt=0, adj=1, xpd=TRUE, cex=0.45)
            par(xpd=TRUE)
            if (length(rnrow) > 50){txtcex = sub.ssdf$cex/2.5 } else { txtcex = sub.ssdf$cex/1.5}
            points(y=sub.ssdf$row, 
                   x=NPLOT - as.numeric(sub.ssdf$node) + 1, 
                   xlim=c(.5,NPLOT+.5), ylim=ylim,
                   ylab='', xlab='', pch=22, cex=txtcex, lwd=.25,
                   col='black', bg=nodecol[nord][as.numeric(sub.ssdf$node)])
            par(xpd=FALSE)
            # Node names:
            xat = 1:NPLOT
            text(y=par()$usr[4] + .005 * diff(par()$usr[3:4]), x=xat, 
                 paste0(1:length(leafrep), '. ', leafrep), col=rev(nodecol[nord]), 
                 srt=90, adj=0, xpd=TRUE, cex=.5)
            box(lwd=.5)
        }

        # pdf(paste0(cadpref, ipref, 'all_go_term.pdf'), width=4, height=6)
        # par(xaxs='i')
        # par(yaxs='i')
        # plot_go_full(1:length(go.ord), leafrep=rev(leafrep[top.nodes[nord]]))
        # dev.off()
    }

    # AGGREGATE ALL (plot small):
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
    print(paste(suid, '----- output with prefix', ipref))
    regfile = paste0(regpref, ipref, '_lreg', suffix)
    load(regfile)
    if (suid == "26920376 - Multiple sclerosis and body mass index (pleiotropy)"){
        ll$log10p = all.regmat[suid,]
        ll$rawlp = all.regmat[suid,]
        ll$df$rawlog10p = all.regmat[suid,]
        ll$df$log10p = all.regmat[suid,]
    }

    if (sum(ll$log10p) > 0){
        finalplot.file = paste0(cadpref, ipref, '_circ_go_snp_panel.pdf')
        if (!file.exists(finalplot.file)){
            print("[STATUS] Plotting page")
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
            # Make the full plot:
            pdf(finalplot.file, width=6, height=8.25)
            par(xaxs='i')
            par(yaxs='i')
            layout(matrix(c(1,2,3,2), ncol=2, byrow=T), widths=c(5,3), heights=c(5, 6), TRUE)
            par(mar=c(1,sp,2,sp))
            circleplot(dend=dd$dend, lab=lab, fractional=FALSE, hit.track=dd$hitdf, 
                       plot.labels=FALSE, add.metadata=TRUE, reduced.metadata=TRUE, nodedf=ntdf, scale=2)
            circos.clear()
            if (s2trait != strait){
                mtext(s2trait, side=3, line=-.5, cex=.5)
                mtext(spmid, side=3, line=-1, cex=.3)
            } else {
                mtext(s2trait, side=3, line=.75, cex=.7)
                mtext(spmid, side=3, line=0, cex=.5)
            }
            mtext(s2init, side=1, line=0, cex=.3)
            smat.reord = smat[snp.cons.rn, rev(nodenames[nord]), drop=F]
            smat.reord = smat.reord[apply(smat.reord,1, sum) > 0,, drop=F]
            plot_snp_full(smat.reord, numsnp=1:nrow(smat.reord), leafrep=rev(leafrep[top.nodes[nord]]))
            if (!(is.null(sub.godf)) && nrow(sub.godf)  > 0){
                plot_go_full(1:length(go.ord), leafrep=rev(leafrep[top.nodes[nord]]))
            }
            dev.off()
        }
    }
}


