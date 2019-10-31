#!/usr/bin/R
# Get list of files corresponding to 127 Roadmap Epigenomes:
source('~/EPIMAP_ANALYSIS/bin/general_EPIMAP_ANALYSIS.R')
library(jsonlite)
library(httr)
today = format(Sys.time(), "%Y%m%d")

# -----------------------
# Get list of accessions:
# -----------------------
rdacc = c()
mainurl = 'https://www.encodeproject.org/'
urlpref = paste0(mainurl, 'report.tsv?type=Experiment&related_series.aliases=roadmap-epigenomics%3A')
eidlist = scan('Annotation/roadmap_eidlist.tsv','c')
for (eid in eidlist){
    eidurl = paste0(urlpref, eid)
    call <- GET(eidurl) 
    df <- read.delim(text=rawToChar(call$content), skip=1)
    rdacc = rbind(rdacc, data.frame(eid = eid, Accession = as.character(df$Accession)))
}

write.table(rdacc, 'Annotation/roadmap_accession_mapping.tsv', quote=F, row.names=F, sep="\t")

# Make a mapping from each Roadmap sample to a BSSID:
accmap = read.delim('accession_bssid_mark.tsv', header=F, stringsAsFactors=F)
names(accmap) = c('Accession','id','Epitope')

mdf = merge(accmap, rdacc, all.x=TRUE)
mdf$eid[is.na(mdf$eid)] = ''

# Ensure all mapped:
missing.eid = eidlist[!(eidlist %in% mdf$eid)]
print(missing.eid)
# E026: Bone Marrow Derived Cultured Mesenchymal Stem Cells (ARCHIVED)
# E049: Mesenchymal Stem Cell Derived Chondrocyte Cultured Cells (ARCHIVED)
arch.expts = rdacc[rdacc$eid %in% missing.eid,]

# GOAL - everything is mapped uniquely.
matching = c()
toresolve = c()
for (eid in eidlist){
    if (!(eid %in% missing.eid)){
        maps = mdf$id[mdf$eid == eid]
        subdf = data.frame(id = maps, eid = eid)
        if (nrow(unique(subdf)) == 1){
            matching = rbind(matching, unique(subdf))
        } else {
            ct = aggregate(eid ~ id, subdf, length)
            names(ct)[2] = 'count'
            ct$eid = eid
            toresolve = rbind(toresolve, ct)
        }
    }
}
# NOTE: Only 75 IDS for some of these tracks.
print(paste("Uniquely matched",nrow(matching), "to", length(unique(matching$id)), 'ids'))

# Resolve mappings:
# 1. Remove all unique mapped already:
eid.resolve = unique(toresolve$eid)
filtdf = toresolve[!(toresolve$id %in% matching$id),]
print(paste(nrow(filtdf),"rows with", length(unique(filtdf$id)), "unique ids"))
# 2. Given unique mappings, take top for each:
filtdf = filtdf[order(-filtdf$count, filtdf$id),]
for (eid in eid.resolve){
    id = filtdf$id[filtdf$eid == eid][1]
    subdf = data.frame(id = id, eid = eid)
    matching = rbind(matching, subdf)
}

# NOTE: Map E026 and E049 to general mesenchymal stem cells (for plots)
matching = rbind(matching, data.frame(id='BSS01261', eid='E026')) # MSC from bone
matching = rbind(matching, data.frame(id='BSS01260', eid='E049')) # MSC from adipose
matching = matching[order(matching$eid),]

write.table(matching, 'Annotation/roadmap_eid_matching.tsv', quote=F, row.names=F, sep="\t")
