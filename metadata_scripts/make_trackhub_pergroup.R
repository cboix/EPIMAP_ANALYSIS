#!/usr/bin/R
# -----------------------------------------------------
# Make the bigwigs trackhubs per group:
# -----------------------------------------------------
domain = system("hostname -d", intern=TRUE)
if (domain == 'broadinstitute.org'){
    bindir='~/data/EPIMAP_ANALYSIS/bin/'
} else {
    bindir='~/EPIMAP_ANALYSIS/bin/'
}
source(paste0(bindir, 'general_EPIMAP_ANALYSIS.R'))
source(paste0(bindir, 'auxiliary_chromImpute_functions.R'))
source(paste0(bindir, 'load_metadata.R'))
library(tidyr)
library(dplyr)
library(stringr)
options(scipen=45) # so we dont get issues writing integers into bedfiles
today <- format(Sys.time(), "%m%d%y")

# -------------------------------
# Load released tracks (metadata)
# -------------------------------
tdf = read.delim('data_to_release/all_released_tracks.tsv', header=F)
names(tdf) = c('id','mark','prefix')
agemap = rd.map[,c('BSSID','Age','AgeUnits')]
names(agemap)[1] = 'id'
meta = merge(meta, agemap, all.x=TRUE)
rownames(meta) = meta$id
meta$fullage = paste0(meta$Age, ' ', meta$AgeUnits, 's')
meta$fullage[grep('unknown',meta$fullage)] = 'unknown'
meta$ctspace = gsub("_", " ", meta$ct)
meta$ctspace = gsub("  ", " ", meta$ctspace)
tdf = merge(tdf, meta, all.x=TRUE)
bwdir = 'https://epigenome.wustl.edu/epimap/data/'
tdf$dir = paste0(bwdir, 'observed/')
tdf$dir[grep('^impute', tdf$prefix)] = paste0(bwdir, 'imputed/')
tdf$dataset = 'Observed'
tdf$dataset[grep('^impute', tdf$prefix)] = 'Imputed'
meta = meta[sort(unique(tdf$id)),]

tdf$mark[tdf$mark == 'DNase-seq'] = 'DNase'
tdf$mark[tdf$mark == 'ATAC-seq'] = 'ATAC'

# -----------------------
# Make subtrack metadata:
# -----------------------
# TODO: Do DNase-seq and Histone marks have diff. properties
tdf$track = paste0("track BW_",  tdf$id, "_", tdf$mark,"_", substr(tdf$dataset,1,3))
tdf$shortLabel = paste("shortLabel", substr(tdf$dataset,1,3), tdf$infoline)
# tdf$longLabel = paste("longLabel", meta$id, meta$ctspace)
tdf$longLabel = paste0("longLabel ", tdf$id, " ", tdf$mark, " ", tdf$name, ' (', tdf$dataset, ')')  # TODO: test if works
tdf$bigDataUrl = paste0('bigDataUrl ', tdf$dir, tdf$prefix, '.bigWig')
tdf$attr = paste(c("type bigWig 0 30", 
                   "viewLimits 0.0:30.0",
                   "noInherit on",
                   "autoScale on",
                   "yLineOnOff on",
                   "maxHeightPixels 64:32:16",
                   "priority 1"), collapse='\n')
tdf$color = sapply(tdf$COLOR, function(x){x = col2rgb(x); paste0('color ',x[1], ',', x[2],',', x[3])})
tdf$metadata = paste('metadata', 
                       paste0('"Epigenome ShortName"="', tdf$infoline,'"'),
                       paste0('"Epigenome FullName"="', tdf$name,'"'),
                       paste0('"Group"="<span style="color:', tdf$COLOR,'">', tdf$GROUP, '</span>"'),
                       paste0('"Type"="', tdf$type,'"'),
                       paste0('"Sex"="', tdf$sex,'"'),
                       paste0('"Age"="', tdf$fullage,'"'),
                       paste0('"Donor"="', tdf$Donor,'"'),
                       sep=" ")
tdf$subGroups = paste0('subGroups view=Sig sampleType=', tdf$id, ' assayType=', tdf$mark, ' dataType=', tdf$dataset, ' donor=', tdf$Donor) 
tdf$parent = paste0('parent EpiMapBWSigView off')

groups = names(colvals$group)
for (group in groups){
    # Making single file:
    groupstr = tolower(gsub("__", "", str_replace_all(group, "[^[:alnum:]]", "_")))
    sdf = tdf[tdf$GROUP == group,]
    track.meta = paste(sdf$track, sdf$shortLabel, sdf$longLabel, sdf$bigDataUrl,
                       sdf$attr, sdf$color, sdf$metadata, sdf$subGroups, sdf$parent, sep="\n")
    print(NROW(sdf))

    genome='hg19'
    hub.meta = c(paste0('hub Epimap_', groupstr),
                 paste0('shortLabel EpiMap release v0.9 (', group, ')'),
                 paste0('longLabel Bigwig tracks for EpiMap v0.9 (', group,')'),
                 'useOneFile on',
                 'email cboix@mit.edu',
                 '\n', 
                 paste('genome', genome))

    # Assays:
    adf = aggregate(id ~ mark, sdf, length)
    adf = adf[order(adf$id, decreasing=T),]
    atypes = adf$mark
    # Donors:
    ddf = aggregate(id ~ Donor, sdf, length)
    ddf = ddf[order(ddf$id, decreasing=T),]
    dtypes = ddf$Donor

    # Overall metadata:
    ovl.meta = c('track EpiMapBWSig',
                 'compositeTrack on',
                 'shortLabel All Signal Tracks',
                 'longLabel All observed and imputed signal tracks from EpiMap',
                 'group epimap',
                 'subGroup1 view Views Sig=Signal',
                 paste(c('subGroup2 sampleType Sample_Type', paste0(meta$id, '=', meta$ct)), collapse=' '),
                 paste(c('subGroup3 assayType Assay_Type ', paste0(atypes, '=', atypes)), collapse=' '),
                 'subGroup4 dataType Data_Type Observed=Observed Imputed=Imputed',
                 paste(c('subGroup5 donor Donor ', paste0(dtypes, '=', dtypes)), collapse=' '),
                 'dimensions dimensionX=assayType dimensionY=sampleType dimA=dataType dimB=donor',
                 'filterComposite on',
                 'sortOrder cellType=+ assayType=+ donor=+',
                 'dragAndDrop on',
                 'visibility hide',
                 'priority 21',
                 'type bed 3',
                 'noInherit on')

    # View track:
    view.meta = c('track EpiMapBWSigView',
                  'shortLabel Signal',
                  'view Sig',
                  'maxHeightPixels 64:32:16',
                  'parent EpiMapBWSig',
                  'type bigWig',
                  'visibility dense')

    # -----------------------
    # Write out the trackHub:
    # -----------------------
    thfile = paste0('trackhubs/pergroup/trackHub_bwtracks_pergroup_', groupstr, '.txt')

    hub.line = paste(hub.meta, collapse='\n')
    ovl.line = paste(ovl.meta, collapse='\n')
    view.line = paste(view.meta, collapse='\n')
    track.line = paste(track.meta, collapse='\n\n')
    writeLines(c(hub.line, ovl.line, view.line, track.line), con=thfile, sep = "\n\n")

    # -----------------------------
    # Make for (new) WUSTL browser:
    # -----------------------------
    # Base string attributes + color (only for legacy):
    sdf$colorstr = sapply(sdf$COLOR, function(x){x = col2rgb(x); paste0('"pr": ',x[1], ', "pg": ', x[2],', "pb": ', x[3])})
    base.info = '{"type": "native_track", "list": [{"name": "refGene" ,"mode":"full"}]}, {"type":"coordinate_override","coord":"chr9,36329955,chr9,37537411"}'
    # Make the per-mark strings:
    sdf$trackline = paste0("{", 
                           paste0('"name": "', substr(sdf$dataset,1,3), " ", sdf$id, " ", sdf$mark,' (', sdf$ct, ')", '),
                           paste0('"url": "', sdf$dir, sdf$prefix, '.bigWig', '", '),
                           '"type": "bigwig", "height": 30, "mode":1, "qtc":{"anglescale":1, "height":30, "summeth":2, "thtype":1, "thmin":0, "thmax":20, "smooth":3,',
                           sdf$colorstr, '}', 
                           ', "metadata": {',
                           paste0('"Assay": "', sdf$mark,'", '),
                           paste0('"Group": "', as.character(sdf$GROUP),'", '),
                           paste0('"Sample": "', sdf$ct,'", '),
                           paste0('"Lifestage": "', sdf$lifestage,'", '),
                           paste0('"Sex": "', sdf$sex,'", '),
                           paste0('"Type": "', sdf$type,'", '),
                           paste0('"DataType": "', sdf$dataset,'"'), "},",
                           paste0('"options": {"color": "', sdf$COLOR, '"}'),
                           "}")
    json.line = paste("[", base.info, ",\n ", paste(sdf$trackline, collapse=",\n "), "]")
    thfile = paste0('trackhubs/pergroup/wustlhub.fulldataset.pergroup.', groupstr, '.json')
    writeLines(json.line, con=thfile)
}

