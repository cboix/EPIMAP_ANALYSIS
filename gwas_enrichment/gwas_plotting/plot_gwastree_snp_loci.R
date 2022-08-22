#!/usr/bin/R
# ------------------------------------------------------
# Script to plot GWAS loci across the 540 enriched GWAS:
# Updated 08/11/20
# TODO: Add rsid
# qsub -cwd -t 1-108 -l h_vmem=20G -l h_rt=1:00:00 -N plt_gwas_loci -j y -b y -V -r y -o $DBDIR/out/ChromImpute "source activate mv_env; $SFTDIR/miniconda3/envs/mv_env/bin/R --slave -f $BINDIR/plot_gwastree_snp_loci.R"
# qsub -cwd -t 2-25 -l h_vmem=20G -l h_rt=1:00:00 -N plt_gwas_loci -j y -b y -V -r y -o $DBDIR/out/ChromImpute "source activate mv_env; $SFTDIR/miniconda3/envs/mv_env/bin/R --slave -f $BINDIR/plot_gwastree_snp_loci.R"
# ------------------------------------------------------
domain = system("hostname -d", intern=TRUE)
if (length(domain) == 0) domain = ''
if (domain == 'broadinstitute.org'){
    bindir='~/data/EPIMAP_ANALYSIS/bin/'
} else {
    bindir='~/EPIMAP_ANALYSIS/bin/'
}
source(paste0(bindir, 'general_EPIMAP_ANALYSIS.R'))
source(paste0(bindir, 'load_metadata.R'))
source(paste0(bindir, 'auxiliary_chromImpute_functions.R'))
library(stringr)
library(scales)
library(GenomicRanges)

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

parpos = function(axis, shift){
    # NOTE: par()$usr is x1, x2, y1, y2
    if (axis == 1) { # X axis 
        par()$usr[1] - shift * diff(par()$usr[1:2])
    } else { # Y axis
        par()$usr[3] - shift * diff(par()$usr[3:4])
    }
}


# Objects we need:
gwcatfile = 'gwascatalog_may03_2019_noquotes.txt'
gwrdafile = sub(".txt", ".Rda", gwcatfile)
gwrsidfile = sub(".txt", "_rsid.txt", gwcatfile)
gwintrdafile = sub(".txt", "_intersections.Rda", gwcatfile)
load(gwrdafile)
load(gwintrdafile)

ddir = 'DHS_Index_WM201902/'
dpref = 'masterlist_DHSs_733samples_WM20180608_all_coords_hg19'
dmlfile = paste0(ddir, dpref, '.core.srt.txt')
dmlnamfile = paste0(ddir, dpref, '_r200_e0_names.core.srt.tsv')
dmlrdafile = paste0(ddir, dpref, '_enh.core.srt.Rda')
load(dmlrdafile)
dmgr = GRanges(enhdf$chr, IRanges(enhdf$start - tol, enhdf$end + tol), name=enhdf$name)

# Tree:
cdllfile = paste0('consensus_object_', usetree, '_062819.Rdata')
load(cdllfile)

ntmeta.rda = paste0('enhancer_tree_metadata.Rda')
load(ntmeta.rda)

gc()


# Prefixes:
extpref = paste0(usetree, '_e', tol, '_')
if (use.adj){ 
    extpref = paste0(extpref, 'adj_')
}
if (use.strict){
    extpref = paste0(extpref, 'p1_')
}
if (use.onecutoff){
    extpref = paste0(extpref, 'onecut_')
}
if (use.adj){ suffix = '_adj1000_10.Rda' } else { suffix = '.Rda' }
if (use.strict){ suffix = '_adj1000_1.Rda' }
if (use.onecutoff){ cutsuff = '_onecut' } else { cutsuff = '' }

# Significant GWAS:
load(paste0(extpref, 'kept_allgwas_ordered', cutsuff, suffix))

# Colors/metadata with "Multiple":
odf = rdcol  # Copy don't change rdcol
odf = rbind(odf, c('Multiple','grey80','Other'))
clvs = c('Tissue', 'Primary Cell', 'Immune', 'Neural', 'Stem-like', 'Other')
odf$category = factor(odf$category, levels=clvs)
# Sort alpha then group:
odf = odf[order(odf$GROUP),]
odf = odf[order(odf$category),]

# Get the leaf reduced information (bag of words, max 3 terms):
NN = nrow(nodetissue)
leafrep = sapply(1:NN, function(x, max.terms=3){
                     blacklist = c('GLAND', 'TISSUE', 'CELL')
                     x = sub("[0-9]*_", "", declist$dec[[x]])
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

# Restrict to just uids we are evaluating:
qdf = qdf[qdf$uid %in% Znam,]
qdf$snploc = gwgr@ranges@start[qdf$queryHits]
dhsmid = (dmgr@ranges@start + dmgr@ranges@width / 2)
qdf$dhsloc = dhsmid[qdf$subjectHits]
qdf$dist = abs(qdf$snploc - qdf$dhsloc)

# ------------------------------------------------
# Get prioritized genes in tissue-specific manner:
# ------------------------------------------------
# Load extra linking objects:
source(paste0(bindir, 'load_linking_objects.R'))
# Get nearest gene per enh:
tssdf = tssdf[!(tssdf$chr %in%  c('chrY', 'chrM')),c('chr','tss','gene')]
tssgr = GRanges(seqnames=tssdf$chr, IRanges(start=tssdf$tss, end=tssdf$tss), name=tssdf$gene)
tssdf = merge(tssdf, genemap)
eout = nearest(dmgr, tssgr)
enhdf$nearest = tssdf$symbol[eout]

# Updated enhancer ids (86):
locgr = with(locdf, GRanges(seqnames=chr, ranges=IRanges(start=start, end=end), name=name))

# Refresh the enhmap (NOTE: linking ind need to be resolved).
enhind = enhdf$cls
enhmap = rep(0, max(enhind))
enhmap[enhind] = 1:length(enhind)

# ---------------------------------------------------------------------
# Load in precomputed objects from extract_gwastree_snp_intersections.R
# ---------------------------------------------------------------------
# Pre-computed group-specific links:
all.group.rda = 'linking_data/all_links_by_group.Rda'
load(all.group.rda)

# The snp x enhancer intersection dataframes:
se.sum.rda = 'linking_data/all_snpxenhancer_intersections.Rda'
snp.sum.rda = 'linking_data/all_snp_summaries_alone.Rda'
enr.sum.rda = 'linking_data/all_enrichment_summaries_alone.Rda'
filt.sum.rda = 'linking_data/all_snp_enrichment_labs_forfilt.Rda'
load(filt.sum.rda)
load(se.sum.rda)
load(snp.sum.rda)
load(enr.sum.rda)

# Pre-calculate the necessary SNPs to intersect for linking:
prelinkrda = 'linking_data/precalc_gwas_SNP_links.Rda'
load(prelinkrda) # Has all.links data.frame

# The tree node links for all enhancers near lead SNPs (table on website):
all.gwlinked.rda = 'linking_data/all_gwas_SNP_links.Rda'
load(all.gwlinked.rda)
print(dim(gw.linksdf))
print("[STATUS] All linking datasets loaded sucessfully.")

# -----------------------------------------------
# Make per-SNP linking plots for a specific GWAS:
# -----------------------------------------------
# Which tasks to run:
CHUNKSIZE = 10
chunk = as.integer(as.numeric(Sys.getenv("SGE_TASK_ID")))
plotrange = ((chunk - 1)* CHUNKSIZE + 1):min(chunk * CHUNKSIZE, length(Znam))
print(paste("[STATUS] Plotting GWAS ids:", paste0(plotrange, collapse=', ')))
plotrange = 11:540
# TODO: Redo 466
plotrange = 467:540
# Where images + linking df go:
gwasdir = 'public_gwas_locus_img/'
for (j in plotrange) {
    # Want to map to closest active, tissue-specific enhancer
    suid = rev(Znam)[j]
    trait = sub(".* - ", "", suid)
    pmid = sub(" - .*", "", suid)
    traitstr = gsub("'","_", gsub(" ", "_", tolower(trait)))
    traitstrnoparen = gsub("<", "lt", gsub(">","gt",
                                           gsub("/", "_", gsub("\\)" ,"", gsub("\\(","",traitstr)))))
    print(paste(j, suid))
    if (suid == "26974007 - Chronic inflammatory diseases (pleiotropy)"){
        fullid = paste0("26974007 - Chronic inflammatory diseases (ankylosing spondylitis, ",
                        "Crohn's disease, psoriasis, primary sclerosing cholangitis, ulcerative colitis) (pleiotropy)")
    } else { fullid = suid }
    # Make UID specific directory for images:
    traitdir = paste0(gwasdir, traitstrnoparen, '_', pmid, '/')
    cmd = paste('mkdir -p', traitdir)
    system(cmd)
    # Get the pre-computed snpdf:
    snpdf = snpenhdf[snpenhdf$uid == suid,]
    snpdf = merge(snpdf, gw.enrdf[gw.enrdf$uid == suid, c('node.rank','node.group')], all.x=TRUE)
    # Initialize data frame for plotted links:
    gwas.plt.linksdf = c()

    # Choose the top groups to add the general-purpose links:
    ndf = unique(snpdf[,c('enrlab','node.group','enr.p')])
    ndf = aggregate(enr.p ~ node.group, ndf, function(x){sum(-log10(x))})
    ndf = ndf[order(ndf$enr.p, decreasing=T),]
    ndf = ndf[!(ndf$node.group %in% c('Multiple', 'Cancer', 'Other')),]
    add.groups = head(ndf$node.group, 10)
    if (nrow(ndf) > 0){
        ldf = ndf[ndf$enr.p >= max(ndf$enr.p) /2 ,]
        link.groups = head(ldf$node.group,4)
    } else { link.groups = c() }
    # NOTE: Pad with diverse groups if we don't have 10 groups:
    add.groups = unique(c(add.groups, 'HSC & B-cell','Muscle','Liver','Brain','Kidney',
                          'Pancreas','Heart', 'Stromal', 'ESC', 'Blood & T-cell'))[1:10]

    # Will make one plot per SNP, get the genes and enhancer names we need:
    uqsnpdf = unique(snpdf[, c('chr','snpPos','snplab')])
    uqsnpdf$loc = paste0(uqsnpdf$chr, '_', uqsnpdf$snpPos)
    uqsnpdf = uqsnpdf[order(uqsnpdf$chr),]
    window = 1e6
    snpgr = with(uqsnpdf, GRanges(seqnames=chr, ranges=IRanges(start=snpPos - window/2, end=snpPos + window/2)))
    uidovldf = as.data.frame(findOverlaps(snpgr, locgr))
    enam = locgr$name[unique(uidovldf$subjectHits)]
    tssovldf = as.data.frame(findOverlaps(snpgr, tssgr))
    egene = tssgr$name[unique(tssovldf$subjectHits)]

    # Get the nearest genes to the GWAS loci:
    sub.gwdf = unique(gwdf[gwdf$uid == fullid, c('chrom','chromStart','pValue')])
    sub.gwdf = sub.gwdf[order(sub.gwdf$chrom),]
    sub.gwdf$chr = paste0('chr', sub.gwdf$chrom)

    # Get the specific links for enhancers around the GWAS lead SNPs:
    sub.linksdf = gw.linksdf[gw.linksdf$uid == suid,]
    sub.linksdf = merge(sub.linksdf, gw.enrdf[gw.enrdf$uid == suid, c('node.rank','node.group')])

    # Get the group-specific links, which we can plot indptly or as background (dashed/thinner lines):
    # TODO: should we subset locdf ?
    sub.enhdf = enhdf[enhdf$name %in% enam,]
    print("[STATUS] Getting group-relevant links:")
    brdf = c()
    for (group in link.groups){
        print(group)
        groupstr = tolower(gsub("__", "", str_replace_all(group, "[^[:alnum:]]", "_")))
        sub.group.rda = paste0('linking_data/links_by_group.', groupstr, '.Rda')
        load(sub.group.rda)
        gdf = merge(gdf, sub.enhdf)
        brdf = rbind(brdf, gdf)
    }
    brdf$GROUP = brdf$group

    # Get all the kept names:

    # Slice locations from hdf5 (matind refers to locdf):
    locdf$mind = 1:nrow(locdf)
    kept.all.locdf = locdf[locdf$name %in% enam,]

    # --------------------------------
    # For each SNP, plot neighborhood:
    # --------------------------------
    tssgr = GRanges(seqnames=tssdf$chr, IRanges(start=tssdf$tss, end=tssdf$tss), name=tssdf$gene)
    center.snps = unique(uqsnpdf$loc)
    cpairsdf = NULL
    for (k in 1:length(center.snps)){
        snp = center.snps[k]
        cat(k, '\t', snp, '\n')
        infoline = uqsnpdf[uqsnpdf$loc == snp,][1,]
        snp.pval = sub(',.*', '', sub('.* \\(', '', infoline$snplab))
        chrom = infoline$chr
        snpgr <- with(infoline, GRanges(seqnames=chr, ranges=IRanges(start=snpPos - window/2, end=snpPos + window/2)))
        gwloc = sub.gwdf$chromStart[sub.gwdf$chr == chrom]

        # Get genes in this locus:
        anno.ovl = data.frame(findOverlaps(snpgr, tssgr))
        kept.gene = as.character(tssgr$name[anno.ovl$subjectHits])
        # kept.gene = kept.gene[kept.gene %in% rownames(gmat)]
        gloc = tssdf$tss[anno.ovl$subjectHits]
        gnames = as.character(sapply(kept.gene, function(x){as.character(genemap$symbol[genemap$gene == x])[1]}))
        ord = order(gloc) 
        gnames = gnames[ord]
        kept.gene = kept.gene[ord]
        gloc = gloc[ord]
        names(gloc) = kept.gene
        # All enhancers in this locus:
        loc.ovl = data.frame(findOverlaps(snpgr, locgr))
        kept.enam = locgr$name[unique(loc.ovl$subjectHits)]

        # Load tested pairs:
        if (is.null(cpairsdf) || cpairsdf$chrom[1] != chrom){
            print(paste("Loading links for", chrom))
            cfile = paste0(lddir, chrom, '_pairs_ENH_ovl_df.tsv.gz')
            cpairsdf = read.delim(cfile, header=T, sep="\t", stringsAsFactors=F)
            cpairsdf$mind = cpairsdf$mind + 1  # Was 0-indexed
            cpairsdf$chrom = chrom
            cpairsdf$name = locdf$name[cpairsdf$mind]

            # Load all for this chromosome:
            kept.chr.locdf = kept.all.locdf[kept.all.locdf$chr == chrom,]
            kept.chr.dhs = kept.chr.locdf$mind
            kept.chr.names = kept.chr.locdf$name

            h5f = H5Fopen(H3K27ac.file)
            h5d = h5f&"matrix"
            kept.mmat = h5d[,matind[kept.chr.dhs]]
            H5Dclose(h5d)
            H5Fclose(h5f)

            # Make average dhs tracks:
            kept.mmat = t(kept.mmat)
            colnames(kept.mmat) = mnames
            tform = make.tform(meta[mnames,'GROUP'], norm=TRUE, u=odf$GROUP)
            kept.avg.mmat = kept.mmat %*% tform 
        }
        # Subset to kept enhancers only:
        sub.ovldf = cpairsdf[cpairsdf$name %in% kept.enam,]

        if (nrow(sub.ovldf) > 0){
            # Get all links:
            sub.ovldf$chr = chrom
            sub.ovldf$mid = locdf$mid[sub.ovldf$mind]
            sub.ovldf = sub.ovldf[order(sub.ovldf$mid),]
            sub.ovldf$dist = sub.ovldf$mid - infoline$snpPos
            sub.ovldf = sub.ovldf[abs(sub.ovldf$dist) < (window / 2),]
            sub.ovldf = sub.ovldf[order(sub.ovldf$mid),]

            # Kept locations:
            kept.locdf = unique(sub.ovldf[,c('mid','mind', 'name')])
            kept.dhs = kept.locdf$mind
            kept.mid = kept.locdf$mid
            kept.sub.names = kept.locdf$name
            rownames(kept.locdf) = kept.locdf$name
            kept.locdf$id = 1:nrow(kept.locdf)

            # Map:
            dind = sapply(kept.sub.names, function(x){which(kept.chr.names == x)})
            sub.avg.mmat = kept.avg.mmat[dind,]

            # Average expression:
            tform = make.tform(meta[colnames(gmat),'GROUP'], norm=TRUE, u=odf$GROUP)
            sub.gmat = gmat[kept.gene,]
            sub.avg.expr = sub.gmat %*% tform
            sub.avg.expr[is.na(sub.avg.expr)] = 0

            use.grouponly = TRUE
            if (use.grouponly){
                if (!is.null(brdf)){
                    preddf = brdf[brdf$name %in% kept.enam | brdf$gene %in% kept.gene,]
                    # preddf = preddf[preddf$score > 0.04,]
                    preddf = preddf[preddf$score > 0.714 / 2,]
                    preddf = merge(preddf, odf)
                } else { preddf = data.frame(score=c()) }
            } else {
                preddf = all.links[all.links$gene %in% kept.gene,]
                preddf = merge(preddf, meta[,c('id','COLOR','GROUP')])
            }

            if (nrow(preddf) > 0){
                npdf = aggregate(score ~ GROUP + COLOR + gene + name, preddf, length)
                npdf = merge(npdf, locdf)
                npdf$mid = (npdf$start + npdf$end) / 2
                npdf = merge(npdf, tssdf)
                npdf$color = as.character(npdf$COLOR)
                npdf$nam = as.character(npdf$GROUP)
                npdf$distToSNP = sapply(npdf$mid, function(x){min(abs(gwloc - x))})
                npdf$on.snp = npdf$distToSNP < 2500
                sublinkdf = npdf
            } else {
                if (ncol(preddf) > 0){
                    sublinkdf = preddf[,c('score','GROUP','COLOR','gene','name')] 
                } else { sublinkdf = preddf }
            }

            # Add the top-links:
            tldf = sub.linksdf[sub.linksdf$name %in% kept.enam, c('name', 'symbol','score','node.group')]
            tldf = tldf[!is.na(tldf$score),]
            scols = c('gene','tss','mid','score', 'color', 'nam', 'on.snp', 'symbol') 
            if (nrow(tldf) > 0){
                tldf = merge(tldf, genemap)
                tldf = merge(tldf, locdf)
                tldf = merge(tldf, tssdf)
                tldf$GROUP = tldf$node.group
                tldf$nam = tldf$node.group
                tldf = merge(tldf,odf)
                tldf$color = tldf$COLOR
                tldf$on.snp = TRUE
                tldf$nam = as.character(tldf$GROUP)
                if (nrow(sublinkdf) > 0){
                    sublinkdf = rbind(sublinkdf[,scols], tldf[,scols])
                } else {
                    sublinkdf = tldf[,scols]
                }
            }

            # Highlight linked genes (only on-snp):
            onsnpdf = sublinkdf[sublinkdf$on.snp,]
            if (nrow(onsnpdf) > 0){
                onsnpdf = aggregate(score ~ color + nam + symbol, onsnpdf, sum)
                onsnpdf = merge(onsnpdf, aggregate(score ~ symbol, onsnpdf, max))
            }

            # Add the plotted links to a full data.frame:
            if (nrow(sublinkdf) > 0){
                subgpdf = sublinkdf[,scols]
                subgpdf$uid = suid
                subgpdf$loc = snp
                gwas.plt.linksdf = rbind(gwas.plt.linksdf, subgpdf)
            }

            plot.red = FALSE
            if (plot.red){
                if (j == 386){
                    # Prostate Cancer:
                    sublinkdf = sublinkdf[sublinkdf$nam %in% c('Cancer','Reproductive'),]
                } else if (suid == '26198764 - Schizophrenia'){
                    # For SCZ:
                    sublinkdf = sublinkdf[sublinkdf$nam %in% c('Brain'),]
                } else {
                    # For coronary artery disease
                    sublinkdf = sublinkdf[sublinkdf$nam %in% c('Liver','Heart','Endothelial'),]
                }
            }

            if (nrow(sublinkdf) > 0){
                # Calc link parameters:
                link.center = with(sublinkdf, (tss + mid) / 2)
                link.radius = abs(link.center - sublinkdf$tss)
            }

            # Set up the data limits:
            xlim = with(infoline, c(snpPos - window /2 , snpPos + window / 2))
            # ylim.atac = c(0, max(sub.avg.mmat, na.rm=T))
            ylim.atac = c(0, 20)

            # For axis:
            step = ((window * 2) / 5)
            mx = (floor(xlim[1] / step) + 1) * step
            atpar = seq(mx, xlim[2], by=step)
            lbls = comma_format()(atpar)

            add.genes = TRUE
            if (add.genes){
                # Default first
                genesuf = '_with_tx'
                genes.height = .5 * 1 + .1
                if (length(gnames) > 0){
                    sub.txdf = aggregate(txid ~ ., txdf[txdf$symbol %in% gnames & txdf$chr == chrom,c('chr','start','end','strand','txid', 'symbol'),], function(x){x[1]})
                    # Adjust the start of the gene to match tss - usually longer transcript, missing, but exists (checked on UCSC browser)
                    sub.txdf = merge(sub.txdf, tssdf[,c('symbol','tss')], all.x=TRUE)
                    pind = which((sub.txdf$strand == '+') & (sub.txdf$start > sub.txdf$tss))
                    nind = which((sub.txdf$strand == '-') & (sub.txdf$end < sub.txdf$tss))
                    if (length(pind) > 0){ sub.txdf$start[pind] = sub.txdf$tss[pind] }
                    if (length(nind) > 0){ sub.txdf$end[nind] = sub.txdf$tss[nind] }
                    sub.txdf$tss = NULL
                    # Add padding for the gene name + link color:
                    if (nrow(onsnpdf) > 0){
                        sub.txdf$is.linked = sub.txdf$symbol %in% onsnpdf$symbol
                    } else { sub.txdf$is.linked = FALSE }
                    sub.txdf$pad = 40000 + 20000 * sub.txdf$is.linked
                    sub.txdf = within(sub.txdf, pad.start <- start - pad * (strand == '+'))
                    sub.txdf = within(sub.txdf, pad.end <- end + pad * (strand == '-'))
                    sub.txdf = sub.txdf[order(sub.txdf$pad.start),]
                    # Figure out how many lines we need:
                    stgr = with(sub.txdf, GRanges(seqnames=chr, IRanges(pad.start, pad.end)))
                    cvg = coverage(stgr)
                    nlines = max(cvg)
                    genes.height = .5 * nlines + .1
                    # Order of genes:
                    tx.ovl = data.frame(findOverlaps(stgr, stgr))
                    tx.ovl = tx.ovl[tx.ovl$queryHits != tx.ovl$subjectHits,]
                    # tx.ovl$qn = sub.txdf$symbol[tx.ovl$queryHits]
                    NGENES = nrow(sub.txdf)
                    sub.txdf$line = rep(0, NGENES)
                    poss.lines = rev(1:nlines)
                    for (i in 1:NGENES){
                        if (i %in% tx.ovl$queryHits){
                            pdf = tx.ovl[tx.ovl$queryHits == i,]
                            comp = sort(pdf[,2])
                            notline = c()
                            for (j in comp){
                                if (i > j){
                                    notline = c(notline, sub.txdf$line[j])
                                }
                            }
                            if (length(notline) > 0){
                                l = poss.lines[!(poss.lines %in% notline)][1]
                            } else { l = nlines }
                            sub.txdf$line[i] = l
                        } else {
                            sub.txdf$line[i] = 1
                        }
                    }
                    sub.exdf = merge(exdf, sub.txdf[,c('txid','symbol', 'line')])
                    sub.uxdf = merge(uxdf, sub.txdf[,c('txid','symbol', 'line')])
                }
            } else { genesuf = ''; genes.height = 0}

            if (plot.red){
                genesuf = paste0(genesuf, '_reduced')
            }
            if (use.grouponly){
                genesuf = paste0(genesuf, '_grouponly')
            }

            NCT = length(add.groups)

            top.arc = TRUE
            plot.width=8
            plot.height=1 + 3 * (NCT / 5) + 1 * add.genes * genes.height / 1.5
            plot.height=2 / 1.2
            imgfilename = paste0(traitdir, 'gene_link_vis_', snp, '_w', window, genesuf, '.png')
            if (!file.exists(imgfilename)){
                png(imgfilename, units='in', res=450, width=plot.width, height=plot.height)
                # arc.sect=3
                arc.sect=2.5
                if (top.arc){
                    lmat = matrix(c(NCT+2 + 1 * add.genes, 1:(NCT + 1 + add.genes * 1)), nrow=NCT + 2 + 1 * add.genes, ncol=1)
                    heights = c(arc.sect, rep(1, NCT + 1))
                } else {
                    lmat = matrix(1:(NCT + 1), nrow=NCT + 1 + 1 * add.genes, ncol=1)
                    heights = c(rep(1, NCT), arc.sect)
                }
                if (add.genes) { heights = c(heights, genes.height) }
                layout(lmat, heights=heights)
                par(xaxs='i')
                sp=0
                par(mar=rep(sp, 4))
                for (i in 1:NCT){ 
                    sgroup = add.groups[i]
                    x = kept.locdf$mid
                    y = c(sub.avg.mmat[,sgroup])
                    ylim = ylim.atac
                    plot(x, y, type='n', axes=F, ylab='', xlab='', ylim=ylim, xlim=xlim)
                    # Add highlights:
                    if (i == NCT){ 
                        axis(1, at=atpar, labels=rep('', length(atpar)), lwd=.25, tck=-0.15)
                        text(x=atpar, y=parpos(2, .55 - .2 * plot.red), labels=paste0(chrom, ":", lbls), xpd=NA, cex=.6)
                    }
                    # Add the highlighted linked locations:
                    tol = 500
                    if (nrow(sublinkdf) > 0){
                        link.atac = sublinkdf$mid
                        rect(xleft=link.atac - tol, xright=link.atac + tol, 
                             ybottom=ylim.atac[1], ytop=ylim.atac[2], col='grey90', border=NA)
                    }
                    # abline(v=gwloc, col='slateblue', lwd=.5, lty='dashed')
                    diffpar = seq(-10, 10, by=2.5) * 1e5 
                    abline(v=infoline$snpPos + diffpar, lwd=.5, col='grey50', lty='dashed')
                    abline(v=gloc, lwd=.5, col='red', lty='dashed')
                    # Plot actual track:
                    rect(xleft=x-tol, ybottom=0, xright=x+tol, ytop=y, col=colvals$group[sgroup], border=NA)
                    groupcex = .8
                    text(x=xlim[1] + 0.028 * diff(xlim), y=mean(ylim.atac), sgroup, 
                         font=2, cex=groupcex, col=colvals$group[sgroup], adj=0)
                    box(lwd=.25)
                    arrowcol='grey25'
                    if (i == 1){
                        if (!(add.genes)) {
                            text(gloc + -1 * diff(xlim) * 1e-3, parpos(2, -.4),
                                 labels=gnames, adj=1, col=arrowcol, font=2, cex=.9)
                        }
                        for (dd in diffpar){
                            tstr = paste(dd / 1e6,'Mb')
                            adj = 1.2 * (dd > 0) - .1
                            if (dd == 0){ tstr=paste0('SNP: ', snp.pval); adj=-.05 }
                            text(x=infoline$snpPos + dd, y = ylim[1] + .75 * (diff(ylim.atac)),
                                 labels=tstr, cex=.5, adj=adj, col='grey30')
                        }
                    }
                }
                # TADs and SNPs on the bottom:
                par(mar=rep(sp, 4))
                plot(0, 1, xlim = xlim, ylim=c(0,1), type='n', axes=F)
                # Plot snps:
                gwind = which(gwloc >= infoline$snpPos - window/2 & gwloc <= infoline$snpPos + window/2)
                if (length(gwind) >  0){
                    kept.gwloc = gwloc[gwind]
                } else {
                    kept.gwloc = gwloc[gwloc >= infoline$snpPos - window]
                }
                rect(xleft=kept.gwloc-tol * 2, xright=kept.gwloc + tol * 2, 
                     ybottom=0.8, ytop=1, border=NA, col='slateblue', lwd=1)
                tcex = 0.575
                # tlab = paste(trait, 'GWAS lead SNPs')
                tlab = paste('GWAS lead SNPs')
                twidth = strwidth(tlab, cex=tcex)
                tloc = xlim[1] + .01 * diff(xlim)
                text(tloc, y=.18, labels=tlab,
                     font=2, cex=tcex, col='slateblue', adj=0)
                min.snp = min(kept.gwloc) # First possible snp.
                max.snp = max(kept.gwloc) # First possible snp.
                # Bottom:
                yloc = c(0.1, 0.75)
                segments(x0=tloc + twidth + .008 *  diff(xlim),
                         x1=max.snp - 0.005 *  diff(xlim),
                         y0=yloc[1],y1=yloc[1], col='slateblue', lwd=.5, xpd=TRUE)
                segments(x0=kept.gwloc - 0.005 *  diff(xlim), x1=kept.gwloc,
                         y0=yloc[1],y1=yloc[2], col='slateblue', lwd=.5, xpd=TRUE)
                segments(x0=kept.gwloc - 0.001 * diff(xlim),
                         x1=kept.gwloc + 0.001 * diff(xlim),
                         y0=yloc[2],y1=yloc[2], col='slateblue', lwd=.5, xpd=TRUE)
                if (add.genes){
                    par(mar=rep(sp, 4))
                    plot(0, 1, xlim = xlim, ylim=c(.5,nlines + .5), type='n', axes=F)
                    # box(lwd=.5)
                    genecol = 'midnightblue'
                    tpad = 1000
                    lpad =0.2 # Up-down padding for exons
                    for (i in 1:nrow(sub.txdf)){
                        genecex = 0.45
                        gene.linked = (sub.txdf$symbol[i] %in% onsnpdf$symbol)
                        if (gene.linked){
                            # Color the genes that are linked:
                            linkcol = onsnpdf[onsnpdf$symbol == sub.txdf$symbol[i],][1,'color']
                            linkgroup = onsnpdf[onsnpdf$symbol == sub.txdf$symbol[i],][1,'nam']
                            txlen = strwidth(sub.txdf$symbol[i], cex=genecex)
                                 # Length of the actual text:
                            with(sub.txdf[i,], rect(xleft=(start - tpad/2) * (strand == '+') + (end + tpad/2) * (strand == '-'),
                                                    xright=(start - 3/2 * tpad - txlen) * (strand == '+') + (end + 3/2 * tpad + txlen) * (strand == '-'),
                                                    ybottom=line - lpad * 2.5, ytop=line + lpad * 2.5, col=linkcol, border=NA))
                        }
                        with(sub.txdf[i,], text(x=(start - tpad) * (strand == '+') + (end + tpad) * (strand == '-'),
                                                y=line, labels=symbol, cex=genecex, xpd=TRUE, adj=(strand == '+') * 1, 
                                                col=ifelse(gene.linked, ifelse(linkgroup %in% c('Multiple','Eye','Lung','Other'),'black','white'),'black')))
                        start = with(sub.txdf[i,], start * (strand == '+') + end * (strand == '-'))
                        end = with(sub.txdf[i,], start * (strand == '-') + end * (strand == '+'))
                        multarrows(x0=start, x1=end,
                                   y0=sub.txdf$line[i], y1=sub.txdf$line[i], 
                                   n_arr = round((sub.txdf$end[i] - sub.txdf$start[i]) / 5000),
                                   length=.02,
                                   lwd=.25, col='midnightblue')
                    }
                    # Add exons:
                    rect(xleft=sub.exdf$start, xright=sub.exdf$end,
                         ybottom=sub.exdf$line - lpad, ytop=sub.exdf$line + lpad, 
                         lwd=.25, col='midnightblue', border=NA)
                    rect(xleft=sub.uxdf$start, xright=sub.uxdf$end,
                         ybottom=sub.uxdf$line - lpad /2, ytop=sub.uxdf$line + lpad /2, 
                         lwd=.25, col='midnightblue', border=NA)
                }
                if (nrow(sublinkdf) > 0){
                    if (top.arc){ 
                        par(mar=rep(sp, 4))
                        plot(0, 1, xlim = xlim, ylim=c(0,1), type='n', axes=F)
                    }
                    asp.ratio = plot.width / (plot.height * (arc.sect / (NCT + add.genes * genes.height + top.arc + arc.sect)))
                    mod.link.radius = abs(link.radius) / cos(pi/4)
                    yvar = mod.link.radius / diff(xlim) * asp.ratio / 1.03
                    # yvar = mod.link.radius / diff(xlim) * asp.ratio / (1.03 + 0.04 * add.genes * (genes.height / 1.5)**2)
                    # if (top.arc){ ang1 = pi/4 }
                    draw.arc(x=link.center, y= 1 * (!top.arc) + (2 * (!top.arc) - 1) * (yvar) * (pi/4), xpd=TRUE, 
                             lwd=(sublinkdf$on.snp == TRUE) * .25 + .25,
                             radius=mod.link.radius, angle1=pi * (!top.arc) + pi/4,angle2=pi * (!top.arc)+ pi - pi/4, col=sublinkdf$color)
                }
                dev.off()
            }
        }
    }

    gwas.plt.rda = paste0(traitdir, traitstrnoparen, '_plotted_links_allsnps_w', window, '.Rda')
    gwas.plt.tsv = paste0(traitdir, traitstrnoparen, '_plotted_links_allsnps_w', window, '.tsv.gz')
    save(gwas.plt.linksdf, file=gwas.plt.rda)
    write.table(gwas.plt.linksdf, file=gzfile(gwas.plt.tsv), quote=F, row.names=F, sep="\t")

}





