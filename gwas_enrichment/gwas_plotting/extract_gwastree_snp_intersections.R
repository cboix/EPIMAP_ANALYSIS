#!/usr/bin/R
# ------------------------------------------------------------------------------
# Script to extract the SNPs and enhancers + links in these loci for the website
# Updated 08/07/20
# -----------------------------------------------------------
domain = system("hostname -d", intern=TRUE)
if (length(domain) == 0) domain = ''
if (domain == 'broadinstitute.org'){
    bindir='~/data/EPIMAP_ANALYSIS/bin/'
} else {
    bindir='~/EPIMAP_ANALYSIS/bin/'
}
source(paste0(bindir, 'general_EPIMAP_ANALYSIS.R'))
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

# Load in + process all of the relevant matrices/datasets:
# commandArgs <- function(trailingOnly=TRUE){
# c(usetree, tol, singlematch, plotting.only, use.adj, use.strict, use.onecutoff) }
# source(paste0(bindir, 'load_statistics_gwastree_enrichments.R'))

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


# ----------------------------------------------
# Reduce gwas summary to make table for website:
# ----------------------------------------------
gwss.rda = 'linking_data/gwascatalog_may03_2019_noquotes_summary.Rda'
if (!file.exists(gwss.rda)){
    sdf = gwssdf[gwssdf$uid %in% Znam,c('pubMedID','trait','sampsize','pValue','initSample','replSample','pubDate', 'uid')]
    names(sdf)[3:4] = c('sampleSize', 'numSNPs')
    rownames(sdf) = NULL
    nrow(unique(sdf))
    gwas.summary.df = sdf
    save(gwas.summary.df, file=gwss.rda)
    # NOTE: that some are doubled because of two replication cohorts:
    # nsdf = aggregate(numSNPs ~ trait + pubMedID, sdf, length)
    # merge(sdf, nsdf[nsdf$numSNPs > 1,c('pubMedID','trait')])
} else {
    load(gwss.rda) 
}



# TODO: for each GWAS, pull down the SNPs, get the enh intersections for all the SNPs
# TODO: Calc distance from enhancer in overlap:
# TODO: use the mapping to the node (any underneath??)

# qdf = qdf[qdf$uid %in% Znam,]
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
tssdf = merge(tssdf, genemap)
tssgr = GRanges(seqnames=tssdf$chr, IRanges(start=tssdf$tss, end=tssdf$tss), name=tssdf$gene)
eout = nearest(dmgr, tssgr)
enhdf$nearest = tssdf$symbol[eout]

# Updated enhancer ids (86) TODO resolve:
locgr = with(locdf, GRanges(seqnames=chr, ranges=IRanges(start=start, end=end), name=name))
en = enhdf$name
ln = locdf$name
# 43 and 39 missing from each:
sum(!(en %in% ln))
sum(!(ln %in% en))

# Refresh the enhmap (NOTE: linking ind need to be resolved).
enhind = enhdf$cls
enhmap = rep(0, max(enhind))
enhmap[enhind] = 1:length(enhind)

# Load actual enrichment matrix:
gtdir = "gwas_tree_analysis/"
regdir = paste0(gtdir, "regressions/")
perdir = paste0(gtdir, "permuted_catalogs/")
epref = paste0(usetree, '_e', tol, '_')
regpref = paste0(regdir, epref)
perpref = paste0(perdir, epref)
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

all.regfile = paste0(regpref, apref, '_logreg_all', cutsuff, suffix)
load(all.regfile)


# -----------------------------------------
# Pre-compute or load group-specific links:
# To be used in plotting around GWAS loci.
# -----------------------------------------
all.group.rda = 'linking_data/all_links_by_group.Rda'
if (!file.exists(all.group.rda)){
    # Enhancer mapping:
    enhdf$loc = with(enhdf, paste(chr, start, end, sep="_"))
    nmapdf = enhdf[,c('loc','name')]
    namemap = enhdf$name
    names(namemap) = enhdf$loc
    # Reduce metadata:
    smeta = meta[meta$id %in% cellorder,]
    groups = as.character(odf$GROUP)
    groups = groups[!(groups %in% c('Cancer','Multiple','Other'))]
    group.linksdf = c()
    for (group in groups){
        groupstr = tolower(gsub("__", "", str_replace_all(group, "[^[:alnum:]]", "_")))
        sub.group.rda = paste0('linking_data/links_by_group.', groupstr, '.Rda')
        if (!file.exists(sub.group.rda)){
            relevant.ids = smeta$id[smeta$GROUP == group]
            cat(group, "\t", length(relevant.ids), "\n")
            # TODO: Merge in as we load data:
            t1 = proc.time()
            gdf = c()
            for (i in 1:length(relevant.ids)){
                id = relevant.ids[i]
                cat(paste0(i, '. ', id,"\t"))
                tab = read.delim(paste0('linking_data/predictions/collated/', id, '_collated_pred.tsv.gz'), header=F)
                names(tab) = c('chr','start','end','gene','score','enh')
                tab$loc = with(tab, paste(chr, start, end, sep="_"))
                cat(nrow(tab),"\t")
                tab$name = namemap[tab$loc]
                if (nrow(tab) > 0){
                    tab$id = id
                    gdf = rbind(gdf, tab[,c('gene','score','name','id')])
                }
                t2 = proc.time() - t1
                cat(round(t2[3],1),"s\n")
            }
            gdf = aggregate(score ~ gene + name + id, gdf, mean) 
            gdf = aggregate(score ~ gene + name, gdf, sum) 
            gdf$score = gdf$score / length(relevant.ids)
            gdf$group = group
            save(gdf, file=sub.group.rda)
        } else { 
            cat('[Loading]', paste0(group,"...\n"))
            load(sub.group.rda) 
            group.linksdf = rbind(group.linksdf, gdf)
        }
    }
    print(dim(group.linksdf))
    save(group.linksdf, file=all.group.rda)
} else {
    load(all.group.rda)
}


# Get the snp x enhancer intersection dataframes:
se.sum.rda = 'linking_data/all_snpxenhancer_intersections.Rda'
snp.sum.rda = 'linking_data/all_snp_summaries_alone.Rda'
enr.sum.rda = 'linking_data/all_enrichment_summaries_alone.Rda'
filt.sum.rda = 'linking_data/all_snp_enrichment_labs_forfilt.Rda'
if (!file.exists(se.sum.rda)){
    snpenhdf = c()
    enhlinkdf = c()
    snplist = list()
    enrlist = list()
    linklist = list()
    for (j in 1:length(Znam)) {
        suid = rev(Znam)[j]
        print(paste(j, suid))
        if (suid == "26974007 - Chronic inflammatory diseases (pleiotropy)"){
            fullid = paste0("26974007 - Chronic inflammatory diseases (ankylosing spondylitis, ",
                        "Crohn's disease, psoriasis, primary sclerosing cholangitis, ulcerative colitis) (pleiotropy)")
        } else { fullid = suid }
        lp = all.regmat[fullid,]
        lind = which(lp > 0)
        lind = lind[order(lp[lind])]
        # Top 25 nodes only:
        lind = tail(lind, 25)
        sqdf = qdf[qdf$uid %in% fullid,c('queryHits','subjectHits','dist')]

        # Assign ALL active enhancer for the node:
        snpdf = ldply(lind, function(i){ 
                        x = cdll$cons[[i]] 
                        x = enhmap[x]; 
                        df = sqdf[sqdf$subjectHits %in% x,]
                        df$node = i
                        return(df) }) 


        
        # Locations and enr/snp pvalues:
        snpdf$loc = paste0(gwdf$chrom[snpdf$queryHits],"_", 
                        gwdf$chromStart[snpdf$queryHits])
        snpdf$p = gwdf$pValue[snpdf$queryHits]
        snpdf$enr.p = 10^(-lp[snpdf$node])
        snpdf = merge(snpdf, data.frame(node=rev(lind),
                                        lab=paste0(1:length(lind), '. ', rev(leafrep[lind]))))
        snpdf$name = dmgr$name[snpdf$subjectHits]
        snpdf$snpPos = as.numeric(sub(".*_", "", snpdf$loc))

        # Save the enhancers in vicinity:
        snpdf = merge(snpdf, enhdf[enhdf$name %in% snpdf$name,], all.x=TRUE)
        snpdf = merge(snpdf, data.frame(node=lind, node.name=leafrep[lind], node.rank=rev(1:length(lind)),
                      node.group=nodetissue$GROUP[lind]), all.x=TRUE)

        # Reorder table:
        snpdf = snpdf[order(snpdf$enr.p),] 
        snpdf = snpdf[order(snpdf$p),] 

        # Choose the top groups to add the general-purpose links:
        ndf = unique(snpdf[,c('node','node.group','enr.p')])
        ndf = aggregate(enr.p ~ node.group, ndf, function(x){sum(-log10(x))})
        ndf = ndf[order(ndf$enr.p, decreasing=T),]
        ndf = ndf[!(ndf$node.group %in% c('Multiple', 'Cancer', 'Other')),]
        ldf = ndf[ndf$enr.p >= max(ndf$enr.p) /2 ,]
        link.groups = head(ldf$node.group,4)

        linklist[[suid]] = link.groups

        # Get the unique names + labs (for filtering):
        snpdf$snplab = with(snpdf, paste0(chr,':', comma_format()(snpPos), ' (p=', sprintf('%0.1e', p),', ', nearest,')'))
        snpdf$enrlab = with(snpdf, paste0(lab, ' (p=', sprintf('%0.1e', enr.p),')'))
        uqsnpdf = unique(snpdf[,c('chr','snpPos','p', 'nearest', 'snplab')])
        uqsnpdf = uqsnpdf[order(uqsnpdf$p),]
        uidsnps = uqsnpdf$snplab
        uqlabdf = unique(snpdf[,c('lab','enr.p', 'enrlab')])
        uqlabdf = uqlabdf[order(uqlabdf$enr.p),]
        uidlabs = uqlabdf$enrlab
        # Add to list:
        snplist[[suid]] = uidsnps
        enrlist[[suid]] = uidlabs

        # Rename and add to full table:
        edf = snpdf[, c('chr','snpPos','start','end','dist','nearest',
                        'node.rank', 'node.name','node.group',
                        'p','enr.p', 'snplab','enrlab')]
        edf$pubMedID = sub(" - .*", "", suid)
        edf$trait = sub(".* - ", "", suid)
        edf$uid = suid
        snpenhdf = rbind(snpenhdf, edf)
    }
    print(dim(snpenhdf))

    # Save the SNPs alone:
    gw.snpdf = unique(snpenhdf[,c('pubMedID','trait','snplab','chr', 'snpPos','nearest','p', 'uid')])
    save(gw.snpdf, file=snp.sum.rda)

    # Save the enrichments alone:
    gw.enrdf = unique(snpenhdf[,c('pubMedID','trait','enrlab','node.rank', 'node.name', 'node.group', 'enr.p', 'uid')])
    save(gw.enrdf, file=enr.sum.rda)

    # Save the SNPs x enhancers x enrichments:
    save(snpenhdf, file=se.sum.rda)

    save(snplist, enrlist, linklist, file=filt.sum.rda)
} else {
    load(filt.sum.rda)
    load(se.sum.rda)
    load(enr.sum.rda)
}


# Options: 1) get group links for overall locus or 2) get the snp-specific links in tree nodes or 3) get snp-specific group links 
# ----------------------------------------------------------
# Pre-calculate the necessary SNPs to intersect for linking:
# ----------------------------------------------------------
prelinkrda = 'linking_data/precalc_gwas_SNP_links.Rda'
if (!file.exists(prelinkrda)){ 
    reqdf = c()
    for (j in 1:length(Znam)) {
        # Want to map to closest active, tissue-specific enhancer
        suid = rev(Znam)[j]
        print(paste(j, suid))
        if (suid == "26974007 - Chronic inflammatory diseases (pleiotropy)"){
            fullid = paste0("26974007 - Chronic inflammatory diseases (ankylosing spondylitis, ",
                            "Crohn's disease, psoriasis, primary sclerosing cholangitis, ulcerative colitis) (pleiotropy)")
        } else { fullid = suid }
        lp = all.regmat[fullid,]
        lind = which(lp > 0)
        lind = lind[order(lp[lind])] # Sort
        lind = tail(lind, 25)
        sqdf = qdf[qdf$uid %in% fullid,c('queryHits','subjectHits','dist')]
        # Assign closest active enhancer for the node:
        snpdf = ldply(lind, function(i){ 
                          x = cdll$cons[[i]] 
                          x = enhmap[x]; 
                          df = sqdf[sqdf$subjectHits %in% x,]
                          # ALL within or just MIN dist:
                          # df2 = aggregate(dist ~ queryHits, df, min)
                          # df = merge(df, df2)
                          df$node = i
                          return(df) }) 
        snpdf$name = enhdf$name[snpdf$subjectHits]
        # Add name x node intersections:
        reqdf = rbind(reqdf, unique(snpdf[,c('node','name')]))
    }

    # Map the indices with the links:
    allids = c()
    nodes = sort(unique(reqdf$node))
    idmap = c()
    for (k in nodes){
        x = declist$dec[[k]] # For breast cancer
        ids = leafmeta[leafmeta$label %in% x,'id']
        idmap = rbind(idmap, data.frame(node=k, id=ids))
    }

    # Load the relevant links:
    linkdf = merge(idmap, reqdf)
    allids = sort(unique(linkdf$id))
    all.links = c()
    t1 = proc.time()
    for (i in 1:length(allids)){
        cat(i,"\t")
        id = allids[i]
        cat(id,"\t")
        tab = read.delim(paste0('linking_data/predictions/collated/', id, '_collated_pred.tsv.gz'), header=F)
        names(tab) = c('chr','start','end','gene','score','enh')
        cat(nrow(tab),"\t")
        # Specific loc:
        enam = linkdf$name[linkdf$id == id]
        sub.enhdf = enhdf[enhdf$name %in% enam,]
        tab = merge(tab, sub.enhdf)
        if (nrow(tab) > 0){
            tab$id = id
            all.links = rbind(all.links, tab[,c('gene','score','nearest', 'enh','name', 'id')])
        }
        t2 = proc.time() - t1
        cat(round(t2[3],1),"s\n")
    }
    save(linkdf, all.links, file=prelinkrda)
} else {
    load(prelinkrda)
}


# Compute the tree node links for all enhancers near lead SNPs:
all.gwlinked.rda = 'linking_data/all_gwas_SNP_links.Rda'
if (!file.exists(all.gwlinked.rda)){
    gw.linksdf = c()
    rownames(enhdf) = enhdf$name
    for (j in 1:length(Znam)) {
        # Want to map to closest active, tissue-specific enhancer
        suid = rev(Znam)[j]
        print(paste(j, suid))
        if (suid == "26974007 - Chronic inflammatory diseases (pleiotropy)"){
            fullid = paste0("26974007 - Chronic inflammatory diseases (ankylosing spondylitis, ",
                            "Crohn's disease, psoriasis, primary sclerosing cholangitis, ulcerative colitis) (pleiotropy)")
        } else { fullid = suid }
        lp = all.regmat[fullid,]
        lind = which(lp > 0)
        lind = lind[order(lp[lind])] # Sort
        lind = tail(lind, 25)
        sqdf = qdf[qdf$uid %in% fullid,c('queryHits','subjectHits','dist')]

        # Assign closest active enhancer for the node:
        snpdf = ldply(lind, function(i){ 
                          x = cdll$cons[[i]] 
                          x = enhmap[x]; 
                          df = sqdf[sqdf$subjectHits %in% x,]
                          df$node = i
                          return(df) }) 
        snpdf$loc = paste0(gwdf$chrom[snpdf$queryHits],"_", 
                           gwdf$chromStart[snpdf$queryHits])
        snpdf$name = enhdf$name[snpdf$subjectHits]

        # Expanded linkdf:
        sldf = unique(merge(linkdf, snpdf))
        nndf = aggregate(id ~ node, unique(sldf[,c('id','node')]), length)
        sldf = unique(merge(all.links, sldf))
        if (nrow(sldf) > 0){
            sldf = aggregate(score ~ id + node + name + gene + dist + subjectHits + loc, sldf, mean)
            sldf = aggregate(score ~ node + name + gene + dist + subjectHits + loc, sldf, sum)
            sldf = merge(sldf, nndf)
            sldf$score = sldf$score / sldf$id
            snpdf = merge(snpdf, sldf, all.x=TRUE)
        } else {
            snpdf = merge(snpdf, sldf[,c('node','score','name','gene','nearest','dist','subjectHits','loc')], all.x=TRUE)
        }
        # Merge back into snpdf, add attributes:
        snpdf = merge(snpdf, enhdf[enhdf$name %in% snpdf$name,], all.x=TRUE)
        snpdf$p = gwdf$pValue[snpdf$queryHits]
        snpdf = merge(snpdf, tssdf, all.x=TRUE)
        snpdf$linkdist = with(snpdf, (start + end) / 2 - tss)
        snpdf$nearest = enhdf[snpdf$name, 'nearest']
        snpdf$enr.p = 10^(-lp[snpdf$node])
        snpdf = merge(snpdf, data.frame(node=lind, node.name=leafrep[lind], node.rank=rev(1:length(lind))))
        # Reorder table:
        snpdf = snpdf[order(snpdf$enr.p, decreasing=T),] 
        snpdf = snpdf[order(snpdf$p),] 
        # Rename and add to full table:
        # TODO: Redo, keeping name of enhancer.
        scols = c('node.rank', 'node.name', 'enr.p', 'loc','p', 'name','dist', 'nearest','symbol','score','linkdist')
        subdf = snpdf[,scols]
        subdf$uid = suid
        gw.linksdf = rbind(gw.linksdf, subdf)
    }
    gw.linksdf$chr = paste0('chr',sub("_.*","",gw.linksdf$loc))
    gw.linksdf$snpPos = as.numeric(sub(".*_","",gw.linksdf$loc))
    gw.linksdf$snpPos = round(gw.linksdf$snpPos)
    gw.linksdf$enr.p = signif(gw.linksdf$enr.p, 2)

    # Reduce, add labs for filtering:
    gw.linksdf$snplab = with(gw.linksdf, paste0(chr,':', comma_format(accuracy=1)(snpPos), ' (p=', sprintf('%0.1e', p),', ', nearest,')'))
    gw.linksdf$enrlab = with(gw.linksdf, paste0(node.rank, '. ', node.name, ' (p=', sprintf('%0.1e', enr.p),')'))

    # Subset and order:
    scols = c('chr','snpPos','p','dist','nearest','symbol','score','linkdist','node.rank','node.name', 'enr.p','name', 'uid', 'enrlab','snplab')
    gw.linksdf = gw.linksdf[,scols]
    gw.linksdf = gw.linksdf[order(gw.linksdf$score),] 
    gw.linksdf = gw.linksdf[order(gw.linksdf$enr.p),] 
    gw.linksdf = gw.linksdf[order(gw.linksdf$p),] 
    save(gw.linksdf, file=all.gwlinked.rda)
} else {
    load(all.gwlinked.rda)
}
print(dim(gw.linksdf))





# ------------------------------------------------------------------------
# Load in the gene expression matrix, get expression for the linked genes:
# ------------------------------------------------------------------------





# -----------------------------------------------
# Make per-SNP linking plots for a specific GWAS:
# -----------------------------------------------
# source(paste0(bindir, 'load_linking_objects.R'))
# enhind = enhdf$cls
# enhmap = rep(0, max(enhind))
# enhmap[enhind] = 1:length(enhind)


j = which(rev(Znam) == '29212778 - Coronary artery disease')
j = which(rev(Znam) == '29892016 - Prostate cancer')
j = which(rev(Znam) == '26198764 - Schizophrenia')
for (j in 11:length(Znam)) {
    # Want to map to closest active, tissue-specific enhancer
    suid = rev(Znam)[j]
    trait = sub(".* - ", "", suid)
    pmid = sub(" - .*", "", suid)
    traitstr = gsub("'","_", gsub(" ", "_", tolower(trait)))
    traitstrnoparen = gsub("/", "_", gsub("\\)" ,"", gsub("\\(","",traitstr)))
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

    # Choose the top groups to add the general-purpose links:
    ndf = unique(snpdf[,c('enrlab','node.group','enr.p')])
    ndf = aggregate(enr.p ~ node.group, ndf, function(x){sum(-log10(x))})
    ndf = ndf[order(ndf$enr.p, decreasing=T),]
    ndf = ndf[!(ndf$node.group %in% c('Multiple', 'Cancer', 'Other')),]
    add.groups = head(ndf$node.group, 10)
    ldf = ndf[ndf$enr.p >= max(ndf$enr.p) /2 ,]
    link.groups = head(ldf$node.group,4)
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
    sub.gwdf = unique(gwdf[gwdf$uid == suid, c('chrom','chromStart','pValue')])
    sub.gwdf = sub.gwdf[order(sub.gwdf$chrom),]
    sub.gwdf$chr = paste0('chr', sub.gwdf$chrom)

    # Get the specific links for enhancers around the GWAS lead SNPs:
    sub.linksdf = gw.linksdf[gw.linksdf$uid == suid,]
    sub.linksdf = merge(sub.linksdf, gw.enrdf[gw.enrdf$uid == suid, c('node.rank','node.group')])

    # Get the group-specific links, which we can plot indptly or as background (dashed/thinner lines):
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

        # Get genes in this locus:
        anno.ovl = data.frame(findOverlaps(snpgr, tssgr))
        kept.gene = as.character(tssgr$name[anno.ovl$subjectHits])
        # kept.gene = kept.gene[kept.gene %in% rownames(gmat)]
        gloc = tssdf$tss[anno.ovl$subjectHits]
        gnames = as.character(sapply(kept.gene, function(x){genemap$symbol[genemap$gene == x]}))
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
            rownames(kept.locdf) = kept.locdf$name
            kept.locdf$id = 1:nrow(kept.locdf)

            # Slice locations from hdf5:
            h5f = H5Fopen(H3K27ac.file)
            h5d = h5f&"matrix"
            sub.mmat = h5d[,matind[kept.dhs]]
            H5Dclose(h5d)
            H5Fclose(h5f)
            sub.mmat = t(sub.mmat)
            colnames(sub.mmat) = mnames
            # Make average dhs tracks:
            tform = make.tform(meta[mnames,'GROUP'], norm=TRUE, u=odf$GROUP)
            sub.avg.mmat = sub.mmat %*% tform 

            # Average expression:
            tform = make.tform(meta[colnames(gmat),'GROUP'], norm=TRUE, u=odf$GROUP)
            sub.gmat = gmat[kept.gene,]
            sub.avg.expr = sub.gmat %*% tform
            sub.avg.expr[is.na(sub.avg.expr)] = 0

            use.grouponly = TRUE
            # TODO: Return plots by two different score thresholds:
            if (use.grouponly){
                preddf = brdf[brdf$name %in% kept.enam | brdf$gene %in% kept.gene,]
                # preddf = preddf[preddf$score > 0.04,]
                preddf = preddf[preddf$score > 0.714 / 2,]
                preddf = merge(preddf, odf)
            } else {
                preddf = all.links[all.links$gene %in% kept.gene,]
                preddf = merge(preddf, meta[,c('id','COLOR','GROUP')])
            }

            # TODO: Reduce properly:
            if (nrow(preddf) > 0){
                npdf = aggregate(score ~ GROUP + COLOR + gene + name, preddf, length)
                npdf = merge(npdf, locdf)
                npdf$mid = (npdf$start + npdf$end) / 2
                npdf = merge(npdf, tssdf)
                npdf$color = as.character(npdf$COLOR)
                npdf$nam = as.character(npdf$GROUP)
                # TODO: TESTING PLOT LIKE THIS:
                npdf$distToSNP = sapply(npdf$mid, function(x){min(abs(gwloc - x))})
                npdf$on.snp = npdf$distToSNP < 2500
                sublinkdf = npdf
            } else {
                sublinkdf = preddf[,c('score','GROUP','COLOR','gene','name')] 
            }

            # Add the top-links:
            tldf = sub.linksdf[sub.linksdf$name %in% kept.enam, c('name', 'symbol','score','node.group')]
            tldf = tldf[!is.na(tldf$score),]
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
                scols = c('gene','tss','mid','score', 'color', 'nam', 'on.snp', 'symbol') 
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
            # GWAS SNPs:
            gwloc = sub.gwdf$chromStart[sub.gwdf$chr == chrom]
            # TODO: Find linked SNPs - for each gene (using gwdf)
            # TODO: FINE MAP each SNP to the closest RELEVANT ENHANCER 

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
                genes.height = .5 * nlines + .1
                if (length(gnames) > 0){
                    sub.txdf = aggregate(txid ~ ., txdf[txdf$symbol %in% gnames,c('chr','start','end','strand','txid', 'symbol'),], function(x){x[1]})
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
                # add.groups = c('Heart','Sm. Muscle','Endothelial','HSC & B-cell')
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
            png(paste0(traitdir, 'gene_link_vis_', snp, '_w', window, genesuf, '_testing.png'), units='in', res=450, width=plot.width, height=plot.height)
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
            # TODO: Add rsid
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
                         # TODO: Switch to highlight direct-to-snp links.
                         lwd=(sublinkdf$on.snp == TRUE) * .25 + .25,
                         radius=mod.link.radius, angle1=pi * (!top.arc) + pi/4,angle2=pi * (!top.arc)+ pi - pi/4, col=sublinkdf$color)
            }
            dev.off()

        }
    }

}





