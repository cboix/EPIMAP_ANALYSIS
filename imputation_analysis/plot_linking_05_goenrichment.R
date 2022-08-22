#!/usr/bin/R
# -----------------------------------
# GO Enrichments to validate linking:
# -----------------------------------
# GO Enrichments linking.
# 1. 
#     a) Get lists of GWAS and their SNPs, 
#     b) Map enhancers to GWAS SNPs
#     c) Highlight the top enriched samples for each GWAS.
# 2. 
#     Then, for each GWAS, you can get a gene list based on which genes the enhancers are linked to (in the top enriched sample). 
#     You can do this for each linking set pretty easily. 
#     You should also do a baseline of linking each enhancer to the nearest gene (that's what GREAT does, a tool for genomic region enrichment analysis)

# 3. Using a gene enrichment tool, like gprofiler (R package: gprofiler2) or similar, you can run GO or other pathway enrichments for the gene set from each linking prediction.

# 4. Starting with one GWAS we understand, we can compare enrichment scores across linking metrics. Then if that makes sense, we can try to branch out to the other GWAS.
# -----------------------------------
domain = system("hostname -d", intern=TRUE)
if (length(domain) == 0) domain = ''
if (domain == 'broadinstitute.org'){
    bindir='~/data/EPIMAP_ANALYSIS/bin/'
} else {
    bindir='~/EPIMAP_ANALYSIS/bin/'
}
source(paste0(bindir, 'general_EPIMAP_ANALYSIS.R'))
source(paste0(bindir, 'auxiliary_chromImpute_functions.R'))
source(paste0(bindir, 'load_metadata.R'))

library(rhdf5)
library(ggplot2)
library(ggpubr)
library(GenomicRanges)
library(Matrix)
library(scales)
library(plotrix)
library(PRROC)
library(huxtable)
library(gprofiler2)

#' Calculate position relative to par()$usr 
#'
#' @param axis 1 is x; 2 is y;
#' @param shift percent of axis to shift over
#' @return position shifted from start of x or y axis
#' @export
parpos = function(axis, shift){
    # NOTE: par()$usr is x1, x2, y1, y2
    if (axis == 1) { # X axis 
        par()$usr[1] - shift * diff(par()$usr[1:2])
    } else { # Y axis
        par()$usr[3] - shift * diff(par()$usr[3:4])
    }
}

# ----------------------------------------
# Make sub-directories for data and plots:
# ----------------------------------------
imgdir = paste0(img, "linking/")
cmd = paste('mkdir -p', imgdir)
system(cmd)
imgpref = paste0(imgdir, 'gwasgo_')

# ---------------------
# Load expression data:
# ---------------------
lddir = 'linking_data/'
tssfile = 'Annotation/gencode.gene.tss.v30lift37.basic.annotation.gtf.gz'
matfile = paste0(lddir, 'merged_log2fpkm.pc.mtx.gz')
generdafile = paste0(lddir, 'gene_mat_tss.Rda')
gmapfile = 'Annotation/gencode.gene.name.mapping.tsv'
if (!file.exists(generdafile)){
    tssdf = read.delim(gzfile(tssfile), header=F, sep="\t", stringsAsFactors=F)
    names(tssdf) = c('chr','tss','tss2','gene')
    tssdf = tssdf[!(tssdf$chr %in%  c('chrY', 'chrM')),]
    # Matrix to wide:
    matlong = read.delim(gzfile(matfile), header=T)
    matwide = spread(matlong, id, log2fpkm)
    gmat = as.matrix(matwide[,-1])
    rownames(gmat) = matwide[,1]
    rm(matwide, matlong)
    # Subset to shared genes:
    kept.genes = sort(tssdf$gene[tssdf$gene %in% rownames(gmat)])
    tssdf = tssdf[tssdf$gene %in% kept.genes,]
    tssdf$gene = factor(tssdf$gene, levels=kept.genes)
    tssdf = tssdf[order(tssdf$gene), ]
    gmat = gmat[kept.genes,]
    gc()
    save(gmat, tssdf, file=generdafile)
} else {
    load(generdafile)
}

genemap = read.delim(gmapfile, header=F)
names(genemap) = c('gene','symbol')
anno = merge(tssdf, genemap) # Note: not same as tssdf file.
annogr = with(anno, GRanges(seqnames=chr, IRanges(tss, tss), name=gene, symbol=symbol))

# --------------------------------------
# Read in the dhs locations and indices:
# --------------------------------------
# TODO: Need the enhancer locations file as well.
lddir = 'linking_data/'
mapfile = 'DHS_Index_WM201902/masterlist_DHSs_733samples_WM20180608_all_coords_hg19_r200_e0_names.core.srt.tsv'
locfile = 'DHS_Index_WM201902/masterlist_DHSs_733samples_WM20180608_all_coords_hg19.core.srt.txt'
indfile = 'masterlist_matindices/matindices_on_mixed_impobs_MODEL_observed_aux_18_on_mixed_impobs_QCUT_intersect_H3K27ac_ENH_masterlist_indices.tsv'
hdfile = 'masterlist_matindices/matindices_on_mixed_impobs_MODEL_observed_aux_18_on_mixed_impobs_QCUT_intersect_H3K27ac_ENH_matrix.hdf5'
namfile = 'masterlist_matindices/matindices_on_mixed_impobs_MODEL_observed_aux_18_on_mixed_impobs_QCUT_intersect_H3K27ac_ENH_matrix_names.tsv'
enames = scan(namfile, 'c') # ENH file is diff. order...
marknamesfile = paste0(lddir, 'mark_matrix_names.txt') # Mark file
dhsrdafile = paste0(lddir, 'dhs_loc.Rda')

# Load data:
if (!file.exists(dhsrdafile)){
    # Locations:
    locmap = read.delim(mapfile, header=T, stringsAsFactors=F)
    locdf = read.delim(locfile, header=F, stringsAsFactors=F)
    colnames(locdf) = c('chr','start','end','name')
    head(locmap)
    head(locdf)
    # Merge mapping:
    locdf = merge(locmap, locdf)
    locdf = locdf[order(locdf$cls),]
    locdf$mid = (locdf$start + locdf$end) / 2
    locmids = locdf$mid  # Keep mids for distances
    rm(locmap)
    gc()
    # Matrix names + indices:
    mnames = scan(marknamesfile,'c')
    matind = as.numeric(scan(indfile, 'c')) + 1
    # Mapping 3.6M to 2M
    enhmap = rep(0,nrow(locdf))
    enhmap[matind] = 1:length(matind)
    locdf = locdf[matind,]
    # Save for loading:
    save(mnames, locdf, matind, enhmap, locmids, file=dhsrdafile)
} else {
    load(dhsrdafile)
}
gc()



# ---------------------------------
# 1a. Assign top sample to GWAS
# ---------------------------------
load('gwas_hg_stats/H3K27ac_merge2_wH3K27ac100_hg_processed_enr.Rda')
cpdf = cplist[['99%']]
# cpdf = calist[['99%']]
colnames(cpdf) = mnames
# sort(scan('gwas_hg_stats/cls_merge2_wH3K27ac100_raw_names.tsv','c'))
cpdf = cpdf[which(apply(cpdf, 1, max) > 0),]
topind = apply(cpdf, 1, which.max)
tdf = data.frame(uid=rownames(cpdf), infoline=meta[colnames(cpdf)[topind], 'infoline'],
                 group=meta[colnames(cpdf)[topind], 'GROUP'], id = colnames(cpdf)[topind])
tdf = tdf[order(tdf$infoline),]
tdf = tdf[order(tdf$group),]
tb = tibble(tdf)
print(tb, n=nrow(tb))

show.mat = FALSE
if (show.mat){
    cmat = cpdf[,which(apply(cpdf, 2, max) > 0)]
    cmat[cmat > 100] = 100
    rmat = reord(cmat)
    rmat = t(reord(t(rmat)))
    rmat = t(diag.mat2(t(rmat))[[1]])
    image(rmat)
}

# ---------------------------------
# 1b. Assign GWAS SNPs to enhancers
# ---------------------------------
# Load the gwas catalog files:
# ----------------------------
gwcatfile = 'gwascatalog_may03_2019_noquotes.txt'
gwrdafile = sub(".txt", ".Rda", gwcatfile)
load(gwrdafile)
# Reduce to just enriched GWAS uids:
kept.uids = unique(tdf$uid)
gwdf = gwdf[gwdf$uid %in% kept.uids,]
gwgr = GRanges(paste0('chr', gwdf$chrom), IRanges(gwdf$chromStart, gwdf$chromEnd), name=gwdf$uid)

gwdf$chr = paste0('chr', gwdf$chrom)

write.table(tdf, 'gwas_mapped_topsample_20200608.tsv', quote=F, row.names=F, sep="\t")
write.table(gwdf[,c('chr','chromStart','chromEnd','uid','trait','pubMedID','pValue')],
            'gwas_mapped_snplocations_20200608.tsv', quote=F, row.names=F, sep="\t")

# --------------------------------------------
# Map kept SNPs to overlapping DHS (ENH only):
# --------------------------------------------
tol = 2500
locgr = with(locdf,  GRanges(chr, IRanges(start - tol, end + tol), cls=cls, name=locdf$name))
ovdf = as.data.frame(findOverlaps(gwgr, locgr))
ovdf$uid = gwgr$name[ovdf$queryHits]
ovdf$cls = locgr$cls[ovdf$subjectHits]
ovdf$name = locgr$name[ovdf$subjectHits]
ovdf$chr = locdf$chr[ovdf$subjectHits]

# Get all tested pairs:
allpairsdf = c()
chrlist = paste0('chr',c(1:22, 'X'))
for (chrom in chrlist){
    print(chrom)
    print(paste("Loading links for", chrom))
    cfile = paste0(lddir, chrom, '_pairs_ENH_ovl_df.tsv.gz')
    cpairsdf = read.delim(cfile, header=T, sep="\t", stringsAsFactors=F)
    cpairsdf$mind = cpairsdf$mind + 1  # Was 0-indexed
    cpairsdf$chrom = chrom
    cpairsdf$cls = matind[cpairsdf$mind]
    cpairsdf = cpairsdf[cpairsdf$cls%in% ovdf$cls ,]
    allpairsdf = rbind(allpairsdf, cpairsdf)
}

# Assign the nearest gene for each GWAS SNP:
nearid = nearest(gwgr, annogr)
neargene = anno[nearid,]

# Load in the H3K27ac corr:
H3K27ac.file = 'linking_data/H3K27ac_all_bin_dense_on_mixed_impobs_r25_e100_allchr_merged.hdf5'
H3K4me1.file = 'linking_data/H3K4me1_all_bin_dense_on_mixed_impobs_r25_e100_allchr_merged.hdf5'
fnames = list.files(path=lddir, pattern='*_precomputed_corr.hdf5')


# For the keratinocyte examples:
# tab = read.delim('linking_data/full_BSS01071_xgpred.txt', header=F, stringsAsFactors=F)
tab = read.delim('linking_data/full_BSS00243_xgpred.txt', header=F, stringsAsFactors=F)
names(tab) = c('chr','start','end','gene','score')
tab = merge(tab, locdf)
tab = aggregate(score ~ cls + gene, tab, max)

# Subset to the top samples --> Could keep all...
# Enhancers: 
suid = "26426971 - Waist-to-hip ratio adjusted for BMI"
suid = '29292387 - Neuroticism (MTAG)'
suid = '23143594 - Psoriasis'
suid = '29212778 - Coronary artery disease'
for (suid in kept.uids){
    print(suid)
    sdf = ovdf[ovdf$uid == suid,]
    # Get any samples with enrichment:
    enr.samp = sort(cpdf[suid, cpdf[suid,] > 0], decreasing=T)
    print(length(enr.samp))
    head(meta[names(enr.samp), c('infoline','GROUP')])

    if (length(enr.samp) > 1){
        sample = names(enr.samp)[1]
    } else {
        sample = names(which(cpdf[suid,] > 0))
    }
    print(meta[sample, c('infoline','GROUP')])
    ematind = which(enames == sample)
    mmatind = which(mnames == sample)

    spairdf = allpairsdf[allpairsdf$cls %in% sdf$cls,]
    spairdf = merge(spairdf, unique(sdf[,c('cls','queryHits')]))
    spairdf = merge(spairdf, locdf)
    spairdf$dist = abs(spairdf$tss - spairdf$mid)
    spairdf$ord = 1:nrow(spairdf)
    cind = spairdf$ind

    # TODO: ONLY WORKS FOR BSS01071 right now.
    stab = merge(spairdf, tab, all.x=TRUE)
    stab$score[is.na(stab$score)] = 0
    stab = stab[order(stab$ord),]

    # Get matrix of correlations:
    # TODO: Compute outside of loop?
    # TODO: Optimize, is quite slow at the moment (may just be too large...)
    mat = matrix(0, ncol=length(fnames), nrow=length(cind))
    marks = c()
    for (i in 1:length(fnames)){
        marks = c(marks, sub("_precomputed.*","",fnames[i]))
        print(marks[length(marks)])
        hfile = paste0(lddir, fnames[i])
        h5f = H5Fopen(hfile)
        mat[, i] = h5f$matrix[cind]
        H5Fclose(h5f)
    }
    colnames(mat) = marks

    h5ls(paste0(lddir, fnames[1]))

    mat[is.na(mat)] = 0 
    # mcoeff = c(.5,1,1,-.5,-.5,.5)
    # cscore = as.numeric(mat %*% as.matrix(mcoeff))

    # H3K27ac values
    h5f = H5Fopen(H3K27ac.file)
    h5d = h5f&"matrix"
    acmark = h5d[mmatind,spairdf$cls]
    H5Dclose(h5d)
    H5Fclose(h5f)
    # H3K4me1 values
    h5f = H5Fopen(H3K4me1.file)
    h5d = h5f&"matrix"
    memark = h5d[mmatind,spairdf$cls]
    H5Dclose(h5d)
    H5Fclose(h5f)

    slist = list()
    slist[['d1']] = spairdf$dist^-1
    slist[['ac_d1']] = acmark * spairdf$dist^-1
    slist[['me_d1']] = memark * spairdf$dist^-1
    slist[['me_corr']] = mat[,'H3K4me1']
    slist[['ac_corr']] = mat[,'H3K27ac'] 
    slist[['me_corr_d1']] = mat[,'H3K4me1']
    slist[['ac_corr_d1']] = mat[,'H3K27ac'] 
    slist[['xg_corr']] = stab$score

    glist = list()
    rlist = list()
    for (score in names(slist)){
        # Score the top gene for each GWAS SNP:
        print(score)
        psdf = data.frame(qh=spairdf$queryHits, gene=spairdf$gene, 
                          val=slist[[score]])
        topdf = aggregate(val ~ qh, psdf, function(x){ if (mean(x) != 0){ max(x) }})
        topdf = merge(topdf, psdf)
        glist[[score]] = topdf$gene

        # Evaluate quality of genes:
        gobj = gost(glist[[score]])
        if (!is.null(gobj$result)){
            rlist[[score]] = gobj$result
            print(paste(score, nrow(rlist[[score]])))
            print(paste(score, max(-log10(rlist[[score]]$p_value))))
        }
    }

    dtab = tibble(rlist[['d1']][,c('term_name','term_size','p_value')])
    rtab = tibble(rlist[['me_corr_d1']][,c('term_name','term_size','p_value')])
    xtab = tibble(rlist[['xg_corr']][,c('term_name','term_size','p_value')])
    print(dtab, n=15)
    print(rtab, n=15)
    print(xtab, n=15)


}





# ---------------------
# Read in the mtx file:
# ---------------------

fullmatfile = 'Enhancer_H3K27ac_matrix_062619.mtx.gz'
enhmatfile = 'Enhancer_H3K27ac_matrix_enhonly_062619.mtx.gz'
enhmargfile = 'Enhancer_H3K27ac_margins_enhonly_062619.tsv.gz'
if (!file.exists(enhmatfile)){
    # Load in the product of enhancer and H3K27ac mtx: 
    mat = read.delim(gzfile(fullmatfile), sep="\t", header=F)
    # Matrix is 0-indexed - turn to 1-indexing
    mat[,1] = mat[,1] + 1
    mat[,2] = mat[,2] + 1




#     b) Map enhancers to GWAS SNPs
#     c) Highlight the top enriched samples for each GWAS. I will send this to you.
# 2. 
#     Then, for each GWAS, you can get a gene list based on which genes the enhancers are linked to (in the top enriched sample). 
#     You can do this for each linking set pretty easily. 
#     You should also do a baseline of linking each enhancer to the nearest gene (that's what GREAT does, a tool for genomic region enrichment analysis)

# ----------------------------
# Load enhancer locations/mat:
# ----------------------------

h5ls(hdfile)
h5f = H5Fopen(hdfile)
h5d = h5f&"matrix"
hcol = h5d[sampind,]
H5Dclose(h5d)
H5Fclose(h5f)


hdfile







# 1a. Assign top sample to GWAS
# 1b. Assign SNPs to GWAS (by sample?)

# load('gwas_hg_stats/H3K27ac_merge2_wH3K27ac100_hg_processed_enr.Rda')
# df = cplist[['99%']]



# -------------------------------------------
# Read and process the GTF for gene plotting:
# -------------------------------------------
parse.gtf.info = function(x){
    x = strsplit(x, "; ")[[1]]
    x = unlist(strsplit(x, " "))
    out = x[seq(2,length(x), by=2)]
    names(out) = x[seq(1,length(x), by=2)] 
    return(out) }

gtf.file = 'Annotation/gencode.v30lift37.basic.annotation.gtf.gz'
gtf.rdafile= 'Annotation/gencode.v30lift37.basic.annotation.Rda'
if (!file.exists(gtf.rdafile)){
    tab = read.delim(gtf.file, sep="\t", skip=5, header=F, stringsAsFactors=F)
    # Find longest tx per gene:
    txdf = tab[tab[,3] == 'transcript',c(1,4,5,7,9)]
    names(txdf) = c('chr','start','end','strand','info')
    tlist = sapply(as.character(txdf$info), parse.gtf.info)
    txdf$txid = sapply(tlist, function(x){x[['transcript_id']]})
    txdf$symbol = sapply(tlist, function(x){x[['gene_name']]})
    txdf$txtype = sapply(tlist, function(x){x[['transcript_type']]})
    txdf$length = txdf$end - txdf$start
    txdf = txdf[txdf$txtype == 'protein_coding',]
    # Reduce to longest transcript only:
    txdf = merge(txdf, aggregate(length ~ symbol, txdf, max))
    txdf$tss = with(txdf,start * (strand == '+') + end * (strand == '-'))
    # Match exons to these transcripts:
    exdf = tab[tab[,3] == 'exon',c(1,4,5,7,9)]
    names(exdf) = c('chr','start','end','strand','info')
    exdf = exdf[grep('protein_coding', exdf$info),]
    elist = sapply(as.character(exdf$info), parse.gtf.info)
    exdf$txid = sapply(elist, function(x){x[['transcript_id']]})
    exdf$exnum = sapply(elist, function(x){x[['exon_number']]})
    exdf$exnum = as.numeric(exdf$exnum)
    exdf = exdf[exdf$txid %in% txdf$txid,]
    exdf = merge(exdf, txdf[,c('txid','symbol')])
    # Match UTRs to these transcripts:
    uxdf = tab[tab[,3] == 'UTR',c(1,4,5,7,9)]
    names(uxdf) = c('chr','start','end','strand','info')
    uxdf = uxdf[grep('protein_coding', uxdf$info),]
    ulist = sapply(as.character(uxdf$info), parse.gtf.info)
    uxdf$txid = sapply(ulist, function(x){x[['transcript_id']]})
    uxdf = uxdf[uxdf$txid %in% txdf$txid,]
    uxdf = merge(uxdf, txdf[,c('txid','symbol')])
    save(exdf, txdf, uxdf, file=gtf.rdafile)
} else {
    load(gtf.rdafile)
}


# -------------------------------------------------------
# Load in validation datasets into genome ranges objects:
# -------------------------------------------------------
vdir = paste0(lddir, 'validation/')
bnames = list.files(path=vdir, pattern='.*Benchmark.*bed')
blist = list()
for (bfile in bnames){
    btype = sub("GM12878.","",sub("-Benchmark.*","", bfile))
    btab = read.delim(paste0(vdir, bfile), header=F)
    names(btab) = c('chr','start','end','gene','score')
    blist[[btype]] = with(btab, GRanges(seqnames=chr, IRanges(start=start, end=end), name=gene, score=score))
    rm(btab)
}

# -------------------------------------
# Read in the precomputed correlations:
# -------------------------------------
# Data files:
H3K27ac.file = 'linking_data/H3K27ac_all_bin_dense_on_mixed_impobs_r25_e100_allchr_merged.hdf5'
H3K4me1.file = 'linking_data/H3K27ac_all_bin_dense_on_mixed_impobs_r25_e100_allchr_merged.hdf5'
fnames = list.files(path=lddir, pattern='*_precomputed_corr.hdf5')

# ------------------------------------------------
# Create some predictions for all chr for BSS00439
# ------------------------------------------------

sample = 'BSS00439'
sampind = which(mnames == sample)
# sampind = 200
meta[mnames[sampind],]
# 1. FDR corr?
# 2. Distance metrics, etc.
# 3. Some mixture

# Get strength for this sample:
# H3K27ac
h5f = H5Fopen(H3K27ac.file)
h5d = h5f&"matrix"
h27samp = h5d[sampind,]
H5Dclose(h5d)
H5Fclose(h5f)

# H3K4me1
h5f = H5Fopen(H3K4me1.file)
h5d = h5f&"matrix"
h4samp = h5d[sampind,]
H5Dclose(h5d)
H5Fclose(h5f)

# Read some xgboost predictions:
# e9tab = read.delim(paste0(lddir, sample, '_E9_2.5.txt'), sep="\t", header=F, stringsAsFactors=F)
# xgdf= read.delim(paste0(lddir, sample, '_E9_0.txt'), sep="\t", header=F, stringsAsFactors=F)
# # xgdf= read.delim(paste0(lddir, sample, '_E9_0.txt'), sep="\t", header=F, stringsAsFactors=F)
# etab = read.delim(paste0(lddir, sample, '_E7_0.txt'), sep="\t", header=F, stringsAsFactors=F)
# xgdf = rbind(xgdf, etab)
# etab = read.delim(paste0(lddir, sample, '_E10_0.txt'), sep="\t", header=F, stringsAsFactors=F)
# xgdf = rbind(xgdf, etab)
xgdf = c()
for (state in c('E7','E8', 'E9','E10', 'E11', 'E15')){
    etab = read.delim(paste0(lddir, sample, '_',state,'_0.txt'), sep="\t", header=F, stringsAsFactors=F)
    xgdf = rbind(xgdf, etab)
}
names(xgdf) = c('chrom','start','end','gene','proba')
xgdf$dhsloc = (xgdf$end + xgdf$start) / 2

chrlist = paste0('chr', c(1:23,'X'))
for (chrom in chrlist){
    # Load links:
    print(paste("Loading links for", chrom))
    cfile = paste0(lddir, chrom, '_pairs_ENH_ovl_df.tsv.gz')
    # genes = sub.gwdf$gene[sub.gwdf$chr == chrom]
    cpairsdf = read.delim(cfile, header=T, sep="\t", stringsAsFactors=F)
    cpairsdf$mind = cpairsdf$mind + 1  # Was 0-indexed
    cpairsdf$chrom = chrom
    cind = cpairsdf$ind

    # Get matrix of correlations:
    mat = matrix(0, ncol=length(fnames), nrow=length(cind))
    marks = c()
    for (i in 1:length(fnames)){
        marks = c(marks, sub("_precomputed.*","",fnames[i]))
        print(marks)
        hfile = paste0(lddir, fnames[i])
        h5f = H5Fopen(hfile)
        mat[, i]= h5f$matrix[cind]
        H5Fclose(h5f)
    }
    colnames(mat) = marks

    mat[is.na(mat)] = 0 
    mcoeff = c(.5,1,1,-.5,-.5,.5)
    cscore = as.numeric(mat %*% as.matrix(mcoeff))

    # Raw values from matrix:
    c27vals = h27samp[matind[cpairsdf$mind]]
    c4vals = h4samp[matind[cpairsdf$mind]]
    # Get strengths, locdf, etc.
    cpairsdf$dhsloc = locdf$mid[cpairsdf$mind]
    cpairsdf$dist = abs(cpairsdf$dhsloc - cpairsdf$tss)


    tol = 100
    tol = 250
    cgr = with(cpairsdf, GRanges(seqnames=chrom, IRanges(start=dhsloc-tol, end=dhsloc+tol), name=gene))
    cpairsdf$acd7 = c27vals * cpairsdf$dist^-.7
    tot.acd7 = aggregate(acd7 ~ gene, cpairsdf, sum)
    names(tot.acd7)[2] = 'tot'
    rownames(tot.acd7) = tot.acd7$gene
    cp2df = cpairsdf[,c('ind','acd7','gene')]
    cp2df$tot = tot.acd7[cp2df$gene, 'tot']
    cp2df$score = cp2df$acd7 / cp2df$tot

    sgdf = xgdf[xgdf$chrom == chrom,]
    sgdf = aggregate(proba ~ dhsloc + gene, sgdf, max)
    sgdf = merge(cpairsdf[,c('gene','ind', 'dhsloc')], sgdf, all.x=TRUE)
    sgdf = sgdf[order(sgdf$ind),]

    ind = which(is.na(sgdf$proba))
    # scfill = cscore * cpairsdf$dist^-1
    # scfill = 0.5 * scfill[ind] / max(scfill[ind])
    # sgdf$proba[ind] = scfill
    sgdf$proba[ind] = 0
    sum(sgdf$proba > 0)

    # Merge xgboost predictions into cpairsdf:
    # SCORES:
    # Distance + distance * strength:
    slist = list()
    slist[['xgboost']] = sgdf$proba
    # slist[['d7']] = cpairsdf$dist^-.7
    slist[['d1']] = cpairsdf$dist^-1
    slist[['xgboost_d1']] = sgdf$proba / (1 - sgdf$proba) * cpairsdf$dist^-1
    # slist[['ac_d7']] = c27vals * cpairsdf$dist^-.7
    slist[['ac_d1']] = c27vals * cpairsdf$dist^-1
    # slist[['me_d7']] = c4vals * cpairsdf$dist^-.7
    # slist[['me_d1']] = c4vals * cpairsdf$dist^-1
    # slist[['all_d7']] = (c27vals + c4vals) * (cpairsdf$dist^-.7)
    # slist[['all2_d7']] = (c27vals + c4vals)^2 * (cpairsdf$dist^-.7)
    # slist[['all_log_d7']] = log10(c27vals + c4vals + 1) * (cpairsdf$dist^-.7)
    # slist[['all_d1']] = (c27vals + c4vals) * (cpairsdf$dist^-1)
    # slist[['cscore_d1']] = cscore * cpairsdf$dist^-1
    # slist[['ac_d7_norm']] = cp2df$score
    slist[['ac_alone']] = c27vals
    # slist[['me_alone']] = c4vals
    slist[['ac_corr']] = mat[,'H3K27ac'] 
    # slist[['me_corr']] = mat[,'H3K4me1']
    # slist[['cscore']] = cscore 
    # slist[['ac_corr_d7']] = mat[,'H3K27ac'] * cpairsdf$dist^-.7
    # slist[['me_corr_d7']] = mat[,'H3K4me1'] * cpairsdf$dist^-.7
    slist[['ac_corr_d1']] = mat[,'H3K27ac'] * cpairsdf$dist^-1
    # slist[['me_corr_d1']] = mat[,'H3K4me1'] * cpairsdf$dist^-1
    # slist[['ac2_corr_d1']] = mat[,'H3K27ac'] * cpairsdf$dist^-1 * c27vals
    # slist[['me2_corr_d1']] = mat[,'H3K4me1'] * cpairsdf$dist^-1 * c4vals

    sapply(slist, length)

    # Colors:
    snap.cols = scan('Annotation/snap_colors.tsv', 'c')

    # MERGE WITH VALIDATIONS
    aucdf = c()
    glist = list()
    for (btype in names(blist)){
        print(btype)
        bgr = blist[[btype]]
        bovl = as.data.frame(findOverlaps(bgr, cgr))
        ind = which(bgr$name[bovl$queryHits] == cgr$name[bovl$subjectHits])
        bovl = bovl[ind,]
        # BG and FG:
        fg = bovl$subjectHits[bgr$score[bovl$queryHits] == 1]
        bg = bovl$subjectHits[bgr$score[bovl$queryHits] == 0]
        print(paste0(length(fg), ', ', length(bg)))

        if (length(fg) > 0){
            curvedf = c()
            for (stype in names(slist)){
                score = slist[[stype]]
                score[is.na(score)] = 0
                roc = roc.curve(scores.class0 = score[fg], scores.class1 = score[bg], curve=TRUE)
                pr = pr.curve(scores.class0 = score[fg], scores.class1 = score[bg], curve=TRUE)
                print(paste0(stype, ': ', round(roc$auc,3), ", ", round(pr$auc.integral,3)))
                aucdf = rbind(aucdf, data.frame(auroc=round(roc$auc,4), auprc=round(pr$auc.integral,3),stype=stype,valid=btype))
                srdf = data.frame(x=roc$curve[,1], y=roc$curve[,2], stype=stype, curve='ROC')
                spdf = data.frame(x=pr$curve[,1], y=pr$curve[,2], stype=stype, curve='PR')
                curvedf = rbind(curvedf, rbind(srdf, spdf))
            }
        }

        gplot = ggplot(curvedf, aes(x,y,col=stype)) + 
            facet_wrap(~curve) + labs(title=btype) + 
            geom_line() + theme_pubr() + 
            theme(legend.position='right')
        ggsave(paste0(imgpref, btype, '_roc_pr_comp_', chrom, '.png'), gplot, units='in', dpi=450, width=9, height=5)
        glist[[btype]] = gplot
    }
    prtab = spread(aucdf[,c('auprc','stype','valid')], valid, auprc)
    roctab = spread(aucdf[,c('auroc','stype','valid')], valid, auroc)

    gplot = patchwork::wrap_plots(glist)
    ggsave(paste0(imgpref, 'allmetric_roc_pr_comp_', chrom, '.png'), gplot, units='in', dpi=450, width=18, height=8)

    prtab = prtab[order(prtab$CHiC, decreasing=T),]
    hp = as_hux(prtab, add_colnames=TRUE)
    position(hp) = "left"
    top_border(hp)[1,]=2
    bottom_border(hp)[1,]=1
    bottom_border(hp)[nrow(hp),]=2
    width(hp) = .75
    quick_pdf(hp, file=paste0(imgpref, "allmetric_auprc_", chrom, ".pdf"))

    roctab = roctab[order(roctab$CHiC, decreasing=T),]
    hp = as_hux(roctab, add_colnames=TRUE)
    position(hp) = "left"
    top_border(hp)[1,]=2
    bottom_border(hp)[1,]=1
    bottom_border(hp)[nrow(hp),]=2
    width(hp) = .75
    quick_pdf(hp, file=paste0(imgpref, "allmetric_auroc_", chrom, ".pdf"))



    # Compare to another validation source (PCHiC)
    # Some validation:
    pcfile = paste0(lddir, 'ActivePromoterEnhancerLinks_JavierrePChIC.tsv.gz')
    pchdf = read.delim(gzfile(pcfile), header=T)
    btdf = unique(pchdf[,c('baitChr','baitSt','baitEnd','baitID')])
    oedf = unique(pchdf[,c('oeChr','oeSt','oeEnd','oeID')])

    # TSS with no window + PCHiC GR objects
    tssgr = GRanges(seqnames=tssdf$chr, IRanges(tssdf$tss, tssdf$tss), name=tssdf$gene)
    btgr = GRanges(seqnames=btdf$baitChr, IRanges(btdf$baitSt, btdf$baitEnd), name=btdf$baitID)
    oegr = GRanges(seqnames=oedf$oeChr, IRanges(oedf$oeSt, oedf$oeEnd), name=oedf$oeID)

    # Overlaps for mapping both names and genes to pchic data:
    btovl = data.frame(findOverlaps(tssgr, btgr))
    btovl$gene = tssgr$name[btovl$queryHits] 
    # Map gene onto pchdf:
    btovl$baitID = btgr$name[btovl$subjectHits] 
    pchdf = merge(pchdf, btovl[,c('gene','baitID')])
    # pchdf$oeMid = (pchdf$oeSt + pchdf$oeEnd) / 2

    # Make all possible links:
    geneuq = unique(pchdf[,c('baitChr','gene', 'baitMid')])
    names(geneuq) = c('chr','gene')
    geneuq = merge(geneuq, anno[,c('gene','tss')])
    genegr = GRanges

    window = 1e6
    gqgr = GRanges(seqnames=geneuq$chr, IRanges(geneuq$tss-window, geneuq$tss+window), name=geneuq$gene)
    govl = data.frame(findOverlaps(gqgr, oegr))
    possdf = data.frame(gene=gqgr$name[govl$queryHits], oeID=oegr$name[govl$subjectHits])
    possdf = merge(possdf, oedf)
    red.pchdf = unique(pchdf[,c('gene','oeID')])
    red.pchdf$hit = 1
    possdf = merge(possdf, red.pchdf, all.x=TRUE)
    possdf$hit[is.na(possdf$hit)] = 0

    btype='PCHiC'
    bgr = with(possdf, GRanges(seqnames=oeChr, IRanges(start=oeSt, end=oeEnd), name=gene, score=hit))
    bovl = as.data.frame(findOverlaps(bgr, cgr))
    ind = which(bgr$name[bovl$queryHits] == cgr$name[bovl$subjectHits])
    bovl = bovl[ind,]
    # BG and FG:
    fg = bovl$subjectHits[bgr$score[bovl$queryHits] == 1]
    bg = bovl$subjectHits[bgr$score[bovl$queryHits] == 0]
    print(paste0(length(fg), ', ', length(bg)))

    if (length(fg) > 0){
        curvedf = c()
        for (stype in names(slist)){
            score = slist[[stype]]
            score[is.na(score)] = 0
            roc = roc.curve(scores.class0 = score[fg], scores.class1 = score[bg], curve=TRUE)
            pr = pr.curve(scores.class0 = score[fg], scores.class1 = score[bg], curve=TRUE)
            print(paste0(stype, ': ', round(roc$auc,3), ", ", round(pr$auc.integral,3)))
            aucdf = rbind(aucdf, data.frame(auroc=round(roc$auc,4), auprc=round(pr$auc.integral,3),stype=stype,valid=btype))
            srdf = data.frame(x=roc$curve[,1], y=roc$curve[,2], stype=stype, curve='ROC')
            spdf = data.frame(x=pr$curve[,1], y=pr$curve[,2], stype=stype, curve='PR')
            curvedf = rbind(curvedf, rbind(srdf, spdf))
        }
    }

    stypes = unique(curvedf$stype)
    scoldf = data.frame(stype=unique(curvedf$stype), 
                        distalone=(stypes %in% c('d1','d7')),
                        ac=(1:length(stypes) %in% grep('ac',stypes)),
                        me=(1:length(stypes) %in% grep('me',stypes)),
                        corr=(1:length(stypes) %in% grep('corr',stypes)),
                        dist=(1:length(stypes) %in% grep('d',stypes)))

    gplot = ggplot(curvedf, aes(x,y,col=stype)) + 
        facet_wrap(~curve) + labs(title=btype) + 
        geom_line() + theme_pubr() + 
        theme(legend.position='right')
    ggsave(paste0(imgpref, btype, '_roc_pr_comp_', chrom, '.png'), gplot, units='in', dpi=450, width=9, height=5)
    glist[[btype]] = gplot


    # -----------------------------------------------------------
    # TEST ALL SAMPLES to check for indexing or validation issue.
    # -----------------------------------------------------------
    # Get strength for this sample:
    # H3K27ac

    btype='CHiC'
    # MERGE WITH VALIDATIONS
    aucdf = c()
    glist = list()
    print(btype)
    bgr = blist[[btype]]
    bovl = as.data.frame(findOverlaps(bgr, cgr))
    ind = which(bgr$name[bovl$queryHits] == cgr$name[bovl$subjectHits])
    bovl = bovl[ind,]
    # BG and FG:
    fg = bovl$subjectHits[bgr$score[bovl$queryHits] == 1]
    bg = bovl$subjectHits[bgr$score[bovl$queryHits] == 0]
    print(paste0(length(fg), ', ', length(bg)))

    for (i in 1:833){
        h5f = H5Fopen(H3K27ac.file)
        h5d = h5f&"matrix"
        ac_mat= h5d[i,]
        H5Dclose(h5d)
        H5Fclose(h5f)
        score = ac_mat[matind[cpairsdf$mind]]
        score[is.na(score)] = 0
        roc = roc.curve(scores.class0 = score[fg], scores.class1 = score[bg], curve=TRUE)
        pr = pr.curve(scores.class0 = score[fg], scores.class1 = score[bg], curve=TRUE)
        print(paste0(i, ': ', round(roc$auc,3), ", ", round(pr$auc.integral,3)))
        aucdf = rbind(aucdf, data.frame(auroc=round(roc$auc,4), auprc=round(pr$auc.integral,3),stype=i,valid=btype))
        srdf = data.frame(x=roc$curve[,1], y=roc$curve[,2], stype=i, curve='ROC')
        spdf = data.frame(x=pr$curve[,1], y=pr$curve[,2], stype=i, curve='PR')
        # curvedf = rbind(curvedf, rbind(srdf, spdf))
    }

    # H3K4me1
    h5f = H5Fopen(H3K4me1.file)
    h5d = h5f&"matrix"
    h4samp = h5d[sampind,]
    H5Dclose(h5d)
    H5Fclose(h5f)



}

