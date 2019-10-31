#!/usr/bin/R
# ---------------------------------------
# Load the distance matrices:
# Imputed, Observed, and Mixed
# Additionally, load in annotation
# ---------------------------------------
domain = system("hostname -d", intern=TRUE)
if (length(domain) == 0) domain = ''
if (domain == 'broadinstitute.org'){
    bindir='~/data/EPIMAP_ANALYSIS/bin/'
} else {
    bindir='~/EPIMAP_ANALYSIS/bin/'
}
library(ComplexHeatmap)
source(paste0(bindir, 'general_EPIMAP_ANALYSIS.R'))
source(paste0(bindir, 'auxiliary_chromImpute_functions.R'))

# Sample to name mapping:
map <- read.delim('Annotation/all_submitted_released_biosample_mapping_20180924.tsv',header=F)
rownames(map) <- map[,2]
names(map) <- c('name','id')

# Track mapping:
tracktab = read.delim('ChromImpute/all_impobs_tracks_uniq_table.tsv', header=F, stringsAsFactors=F)
names(tracktab) = c('uqsample','Epitope','prefix')

# Color mapping:
rdcol = read.delim('Annotation/updated_tissue_colors.tsv', header=T, stringsAsFactors=F)
names(rdcol) = c('GROUP','COLOR', 'category')
rdcol = rdcol[order(rdcol$COLOR, decreasing=T), ]

# Metadata for imputation datasets:
rd.map = read.delim('Annotation/updated_imputation_metadata.tsv',header=T, stringsAsFactors=F)
meta = rd.map[,c('BSSID','sample.name','newgroup','secondary','Extended.Info','origin','perturb','lifestage','age','sex','type')]
names(meta)[1:5] = c('id','ct','GROUP','SECONDARY','infoline')
rdcol = rdcol[rdcol$GROUP %in% unique(meta$GROUP),]
meta = merge(meta, rdcol)
meta[is.na(meta)] <- ''
# NOTE: Only include if we want to add the 10 modified samples:
# rdcol = rbind(rdcol, data.frame(GROUP='Modified',COLOR='lightgrey'))

# Colors:
odf = rdcol  # Copy don't change rdcol
# odf = rbind(odf, c('Multiple','grey80','Other'))
clvs = c('Tissue', 'Primary Cell', 'Immune', 'Neural', 'Stem-like', 'Other')
odf$category = factor(odf$category, levels=clvs)
# Sort alpha then group:
odf = odf[order(odf$GROUP),]
odf = odf[order(odf$category),]

colvals = list()
colvals[['sex']] = c('male' = 'lightblue', 'female' = 'pink', 'unknown/mixed' = 'white')
colvals[['type']] = c('tissue' = 'darkgoldenrod1', 'cell line' = 'royalblue', 'primary cell' = 'indianred', 'in vitro differentiated cells' = 'mediumpurple1', 'unknown/mixed'='white')
colvals[['lifestage']] = c('embryonic' = 'forestgreen', 'adult' = 'royalblue', 'newborn' = 'indianred', 'child' = 'mediumpurple3', 'unknown/mixed' = 'white')
colvals[['project']] = c('ENCODE (New)' = 'darkblue', 'ENCODE 2012'='darkgrey', 'Roadmap 2015' = 'lightgrey', 'Roadmap (New)'='indianred', 'GGR'='forestgreen')
colvals[['group']] = odf$COLOR
names(colvals[['group']]) = odf$GROUP

# PCH:
pchvals = list()
pchvals[['sex']] = c('male'=1, 'female'=19, 'unknown/mixed'=18)
pchvals[['type']] = c('tissue' = 19, 'cell line'=1, 'primary cell'=1, 'in vitro differentiated cells'=1, 'unknown/mixed'=1)
pchvals[['lifestage']] = c('embryonic'=1, 'adult'=19, 'newborn'=19, 'child'=19, 'unknown/mixed'=19)
pchvals[['Project']] = c('ENCODE (New)'=19, 'ENCODE 2012'=1, 'Roadmap 2015'=1, 'Roadmap (New)'=19, 'GGR'=19)

# Get availability from sample/mark table:
smat = read.delim('ChromImpute/sample_mark_table.tsv', header=F)
names(smat) <- c('id','Epitope','file')
smat$file = 1
cells <- aggregate(file ~ id, smat, length)
epitopes <- aggregate(file ~ Epitope, smat, length)

# Matrix of mark availability
wm <- spread(smat, id, file, fill=0)
rownames(wm) <- wm$Epitope 
avail = as.matrix(wm[,-1])
rownames(avail) <- wm$Epitope 

# REMOVE: Cell-types that only have ATAC-seq:
ared = avail[,avail['ATAC-seq',] > 0]
atacIDs = names(which(apply(ared,2,sum) == 1))
dim(tracktab[tracktab$uqsample %in% paste0(atacIDs, '_obs'),])
dim(tracktab[tracktab$uqsample %in% paste0(atacIDs, '_imp'),])


# Availability meta:
eporder = epitopes$Epitope[order(epitopes$file, decreasing=T)]
eporder = as.character(eporder)
cellorder = as.character(cells$id[order(cells$file, decreasing=T)])
main.marks = as.character(head(eporder, 9))
NMAIN = length(main.marks)

# Get availability from sample/mark table:
extmat = read.delim('ChromImpute/sample_mark_table_extended.tsv', header=F)
names(extmat) <- c('id','Epitope','file')
extmat$file=1
# Matrix of mark availability
ewm <- spread(extmat, id, file, fill=0)
rownames(ewm) <- ewm$Epitope 
extavail = as.matrix(ewm[,-1])
rownames(extavail) <- ewm$Epitope 
extepitopes <- aggregate(file ~ Epitope, extmat, length)
exteporder = extepitopes$Epitope[order(extepitopes$file, decreasing=T)]
exteporder = c(as.character(eporder[1:13]), 
               as.character(exteporder[!(exteporder %in% eporder)]), 
               as.character(eporder[14:35]))

# -----------------------
ct.list = cellorder
N <- length(ct.list)
# -------------------------------------------
# Label samples by project of origin:
# Separate samples by: 
# 1. ENCODE 2012 (ENCODE + ROADMAP?)
# 2. Within Roadmap 111 genomes (any overlap)
# 3. ENCODE new
# 4. Roadmap new
# 5. GGR
# -------------------------------------------
# Roadmap annotations:
rd.matching = read.delim('Annotation/roadmap_eid_matching.tsv', stringsAsFactors=F)
roadmap.colors = read.delim('Annotation/roadmap_eid_colormap.tsv', stringsAsFactors=F)
rdacc = read.delim('Annotation/roadmap_accession_mapping.tsv', stringsAsFactors=F)

# Load information on whether a sample is roadmap or not:
accfile = 'accession_bssid_mark_extended.tsv'
if (!file.exists(accfile)){
    accmap = read.delim('accession_bssid_mark.tsv', header=F, stringsAsFactors=F)
    names(accmap) = c('Accession','id','Epitope')
    filemeta = read.delim('Annotation/file_extended_metadata_20180128.tsv', header=T, stringsAsFactors=F)
    accmap = merge(accmap, unique(filemeta[,c('Accession','Project',
                                              'replicates.library.biosample.donor.accession')]), all.x=T)
    names(accmap)[5] = 'Donor'
    accmap$is.ENC = accmap$Project == "ENCODE"
    accmap$is.RDM = accmap$Project == "Roadmap"
    accmap$is.GGR = accmap$Project == "GGR"
    # See if was used in Roadmap 2015:
    accmap$in.RM = sapply(accmap$Accession, function(x){x %in% rdacc$Accession})
    accmap$total = 1
    write.table(accmap, accfile, quote=F, sep="\t", row.names=F)
} else {
    accmap = read.table(accfile, sep="\t", stringsAsFactors=F, header=T)
}

aw = aggregate(cbind(in.RM, is.ENC, is.RDM, is.GGR, total) ~ id , accmap, sum)

# Make labels for data:
aw$Project = 'GGR'
aw$Project[aw$is.ENC > 0] = 'ENCODE (New)'
aw$Project[aw$is.RDM > 0] = 'Roadmap (New)'
aw$Project[aw$in.RM > 0 & aw$is.ENC > 0] = 'ENCODE 2012'
aw$Project[aw$in.RM > 0 & aw$is.RDM > 0] = 'Roadmap 2015'

# Final project labels:
projdf = unique(aw[, c('id','Project')])
meta = merge(meta, projdf, all.x=TRUE)
meta$Project[is.na(meta$Project)] = 'ENCODE (New)'
projecttotals = aggregate(id ~ Project, meta, length)

# Add Donor metadata:
donordf = aggregate(Donor ~ id, accmap, function(x){
                        x = as.character(unlist(sapply(x, function(x){strsplit(x,",")[[1]]})))
                        x = sort(unique(x))
                        paste0(x, collapse=',')})
meta = merge(meta, donordf, all.x=TRUE)
dlist = aggregate(id ~ Donor, meta, length)
dlist = dlist[dlist$id > 1,]
# NOTE: 126 repeated donors

# In case some aren't labeled:
meta = merge(meta, map[ct.list,], all.y=TRUE, all.x=TRUE)
meta$GROUP[is.na(meta$GROUP)] = 'Modified'
meta$COLOR[meta$GROUP == 'Modified'] = 'lightgrey'
# Remove modified:
keep.nomod = which(meta$GROUP != 'Modified')
length(which(meta$GROUP == 'Modified'))
meta = meta[keep.nomod,]
meta$GROUP = factor(meta$GROUP, levels = as.character(rdcol$GROUP))
rownames(meta) = meta$id
meta$type[meta$type == 'NA'] = 'unknown/mixed'
meta$type[is.na(meta$type)] = 'unknown/mixed'
meta$sex[is.na(meta$sex)] = 'unknown/mixed'
meta$sex[meta$sex == 'NA'] = 'unknown/mixed'
meta$sex[!c(meta$sex %in% c('male', 'female'))] = 'unknown/mixed'
meta$lifestage[meta$lifestage == 'NA'] = 'unknown/mixed'
meta$lifestage[!(meta$lifestage %in% c('embryonic', 'adult', 'newborn', 'child'))] = 'unknown/mixed'
meta$lifestage[is.na(meta$lifestage)] = 'unknown/mixed'

# Reduced cell set (no modifieds or atac-only):
remove_poorq = c('BSS00456', 'BSS00457', 'BSS00474', 'BSS01880')
cellorder = cellorder[cellorder %in% meta$id]
cellorder = cellorder[!(cellorder %in% atacIDs)]
cellorder = cellorder[!(cellorder %in% remove_poorq)]

# Metadata matrix for plotting:
metamat = cbind('lifestage'=meta$lifestage, 'sex'=meta$sex,
                'type'=meta$type, 'project'=meta$Project,
                'group'=as.character(meta$GROUP))
meta$type[is.na(meta$type)] = 'ENCODE'
rownames(metamat) = meta$id

palette = colryb

# Mark definition:
mainmarks = c('H3K27ac','H3K27me3','H3K36me3','H3K4me1','H3K4me3','H3K9me3')
tier2marks = c("ATAC-seq", "DNase-seq","H2AFZ", "H3K4me2", "H3K79me2", "H3K9ac", "H4K20me1")
NMAIN = length(mainmarks)
t12marks = c(mainmarks, tier2marks)
NT12 = length(t12marks)

# Punctate vs. broad definitions:
pm = c('ATAC-seq','DNase-seq','H2AFZ','H3K27ac','H3K4me1','H3K4me2','H3K4me3','H3K9ac')
bm = c('H3K27me3','H3K36me3','H3K79me2','H3K9me3','H4K20me1')
markdef = rbind(data.frame(mark=pm, type='punctate'), data.frame(mark=bm, type='broad'))

# For legend:
col.list = list()
col.list.h = list()
for (attr in names(colvals)){
    acol = metamat[cellorder, attr]
    nv = sapply(names(colvals[[attr]]), function(x, mat=acol){ paste0(x, ' (', sum(mat == x), ')') })
    col.list[[attr]] = Legend(labels = nv, legend_gp = gpar(fill = colvals[[attr]]), title = capitalize(attr))
    col.list.h[[attr]] = Legend(labels = nv, legend_gp = gpar(fill = colvals[[attr]]), title = capitalize(attr), direction='horizontal')
}
pd.legend = packLegend(col.list$group, col.list.h$project, col.list$sex, col.list$type, col.list$lifestage)
pd.legend.horiz = packLegend(col.list.h$group, col.list.h$project,col.list.h$sex, col.list.h$type, col.list.h$lifestage, direction='horizontal')

