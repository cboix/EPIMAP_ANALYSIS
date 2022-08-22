#!/usr/bin/R
# -----------------------------------------------------
# Make the bigwigs trackhub:
# Have observed, imputed, and by sample slicing
# Color tracks by either sample type or by obs/imputed?
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
options(scipen=45) # so we dont get issues writing integers into bedfiles
today <- format(Sys.time(), "%m%d%y")


# -------------------------------
# Load released tracks (metadata)
# -------------------------------
tdf = read.delim('data_to_release/all_released_tracks.tsv', header=F)
names(tdf) = c('id','mark','prefix')
agemap = rd.map[,c('BSSID','age','ageUnits')]
names(agemap)[1] = 'id'
meta = merge(meta, agemap, all.x=TRUE)
rownames(meta) = meta$id
meta$fullage = paste0(meta$age, ' ', meta$ageUnits, 's')
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
track.meta = paste(tdf$track, tdf$shortLabel, tdf$longLabel, tdf$bigDataUrl,
                   tdf$attr, tdf$color, tdf$metadata, tdf$subGroups, tdf$parent, sep="\n")

# Assays and donors:
adf = aggregate(id ~ mark, tdf, length)
adf = adf[order(adf$id, decreasing=T),]
atypes = adf$mark

ddf = aggregate(id ~ Donor, tdf, length)
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


# --------------------------------
# Write out the chromHMM trackHub:
# --------------------------------
thfile = 'trackHub_bwtracks_v1.txt'

ovl.line = paste(ovl.meta, collapse='\n')
view.line = paste(view.meta, collapse='\n')
track.line = paste(track.meta, collapse='\n\n')

writeLines(c(ovl.line, view.line, track.line), con=thfile, sep = "\n\n")


# -------------------------------
# Additional -  ByAssay BySample:
# -------------------------------
# By Assay:
super.meta = c('track CompEpiMapbyAssay',
               'superTrack on',
               'shortLabel By Assay',
               'longLabel EpiMap data by assay',
               'group epimap',
               'visibility hide',
               'priority 11')

super.line = paste(c("\n", super.meta, "\n"), collapse='\n')
write(super.line, file=thfile, append=T)

for (assay in atypes){
    print(assay)
    astr = sub("-", "", assay)
    sdf = tdf[tdf$mark == assay,]
    smeta = unique(sdf[,c('id','ct')])
    ddf = aggregate(id ~ Donor, sdf, length)
    ddf = ddf[order(ddf$id, decreasing=T),]
    dtypes = ddf$Donor

    # Assay and assay view lines:
    assay.meta = c(paste0('track CompEpiMap', astr),
                   'compositeTrack on',
                   'parent CompEpiMapbyAssay',
                   paste('shortLabel', assay),
                   paste('longLabel',assay,'tracks for', nrow(smeta), 'sample type(s)'),
                   'group epimap',
                   'subGroup1 view Views ALN=Alignments COV=Coverage',
                   paste(c('subGroup2 sampleType Sample_Type', paste0(smeta$id, '=', smeta$ct)), collapse=' '),
                   'subGroup3 dataType Data_Type Observed=Observed Imputed=Imputed',
                   paste(c('subGroup4 donor Donor ', paste0(dtypes, '=', dtypes)), collapse=' '),
                   'dimensions dimensionX=dataType dimensionY=sampleType',
                   'sortOrder sampleType=+ dataType=+',
                   'dividers sampleType',
                   'dragAndDrop on',
                   'priority 1',
                   'visibility dense',
                   'type bed 3')
    assayview.meta = c(paste0('track CompEpiMap', astr, 'ViewCOV'),
                       'shortLabel Coverage',
                       'view COV',
                       'maxHeightPixels 64:32:16',
                       paste0('parent CompEpiMap', astr),
                       'visibility dense')

    assay.line = paste(c("\n", assay.meta, "\n"), collapse='\n')
    assayview.line = paste(c("\n", assayview.meta, "\n"), collapse='\n')
    write(assay.line, file=thfile, append=T)
    write(assayview.line, file=thfile, append=T)

    # -----------------------
    # Make subtrack metadata:
    # -----------------------
    # TODO: Do DNase-seq and Histone marks have diff. properties
    # tdf$track = paste0("track BW_",  tdf$id, "_", tdf$mark,"_", substr(tdf$dataset,1,3))
    sdf$track = paste0("track CompAssayBW_",  sdf$id, "_", sdf$mark,"_", substr(sdf$dataset,1,3))
    sdf$subGroups = paste0('subGroups view=COV sampleType=', sdf$id, ' dataType=', sdf$dataset, ' donor=', sdf$Donor) 
    sdf$parent = paste0('parent CompEpiMap', astr, 'ViewCOV off')
    track.meta = paste(sdf$track, sdf$shortLabel, sdf$longLabel, sdf$bigDataUrl,
                       sdf$attr, sdf$color, sdf$metadata, sdf$subGroups, sdf$parent, sep="\n")


    track.line = paste(track.meta, collapse='\n\n')
    write(track.line, file=thfile, append=T)

}


# ----------
# By Sample:
# ----------
super.meta = c('track CompEpiMapbySample',
               'superTrack on',
               'shortLabel By Sample',
               'longLabel EpiMap data by sample',
               'group epimap',
               'visibility hide',
               'priority 11')

super.line = paste(c("\n", super.meta, "\n"), collapse='\n')
write(super.line, file=thfile, append=T)

samples = unique(tdf$id)

for (sample in samples){
    print(sample)
    sdf = tdf[tdf$id == sample,]
    smeta = unique(sdf[,c('id','ct', 'infoline', 'Donor')])
    ddf = aggregate(id ~ Donor, sdf, length)
    ddf = ddf[order(ddf$id, decreasing=T),]
    dtypes = ddf$Donor

    # Assays and donors:
    adf = aggregate(id ~ mark, sdf, length)
    subassays = atypes[atypes %in% adf$mark]

    # Assay and assay view lines:
    assay.meta = c(paste0('track CompEpiMapSamp', sample),
                   'compositeTrack on',
                   'parent CompEpiMapbySample',
                   paste('shortLabel', sample, smeta$infoline),
                   paste0('longLabel ',sample,' (',smeta$ct, ') tracks for ', length(subassays), ' assay(s)'),
                   'group epimap',
                   'subGroup1 view Views ALN=Alignments COV=Coverage',
                   paste(c('subGroup2 assayType Assay_Type ', paste0(subassays, '=', subassays)), collapse=' '),
                   'subGroup3 dataType Data_Type Observed=Observed Imputed=Imputed',
                   paste(c('subGroup4 donor Donor ', paste0(dtypes, '=', dtypes)), collapse=' '),
                   'dimensions dimensionX=dataType dimensionY=assayType',
                   'sortOrder assayType=+ dataType=+',
                   'dividers assayType',
                   'dragAndDrop on',
                   'priority 1',
                   'visibility dense',
                   'type bed 3')

    assayview.meta = c(paste0('track CompEpiMapSamp', sample, 'ViewCOV'),
                       'shortLabel Coverage',
                       'view COV',
                       'maxHeightPixels 64:32:16',
                       paste0('parent CompEpiMapSamp', sample),
                       'visibility dense')

    assay.line = paste(c("\n", assay.meta, "\n"), collapse='\n')
    assayview.line = paste(c("\n", assayview.meta, "\n"), collapse='\n')
    write(assay.line, file=thfile, append=T)
    write(assayview.line, file=thfile, append=T)

    # -----------------------
    # Make subtrack metadata:
    # -----------------------
    # TODO: Do DNase-seq and Histone marks have diff. properties
    # tdf$track = paste0("track BW_",  tdf$id, "_", tdf$mark,"_", substr(tdf$dataset,1,3))
    sdf$track = paste0("track CompSampBW_",  sdf$id, "_", sdf$mark,"_", substr(sdf$dataset,1,3))
    sdf$subGroups = paste0('subGroups view=COV assayType=', sdf$mark, ' dataType=', sdf$dataset, ' donor=', sdf$Donor) 
    sdf$parent = paste0('parent CompEpiMapSamp', sample, 'ViewCOV off')
    track.meta = paste(sdf$track, sdf$shortLabel, sdf$longLabel, sdf$bigDataUrl,
                       sdf$attr, sdf$color, sdf$metadata, sdf$subGroups, sdf$parent, sep="\n")

    track.line = paste(track.meta, collapse='\n\n')
    write(track.line, file=thfile, append=T)

}



# ----------
# By group:
# ----------
super.meta = c('track CompEpiMapbyGroup',
               'superTrack on',
               'shortLabel By Group',
               'longLabel EpiMap data by group',
               'group epimap',
               'visibility hide',
               'priority 11')

super.line = paste(c("\n", super.meta, "\n"), collapse='\n')
write(super.line, file=thfile, append=T)

groups = sort(as.character(unique(tdf$GROUP)))

for (group in groups){
    print(group)
    groupstr = gsub(" & ","and",group)
    groupstr = gsub("[. -]","",groupstr)
    sdf = tdf[tdf$GROUP == group,]
    smeta = unique(sdf[,c('id','ct', 'infoline', 'Donor')])
    ddf = aggregate(id ~ Donor, sdf, length)
    ddf = ddf[order(ddf$id, decreasing=T),]
    dtypes = ddf$Donor

    # Assays and donors:
    adf = aggregate(id ~ mark, sdf, length)
    subassays = atypes[atypes %in% adf$mark]

    # Assay and assay view lines:
    assay.meta = c(paste0('track CompEpiMapGroup', groupstr),
                   'compositeTrack on',
                   'parent CompEpiMapbyGroup',
                   paste('shortLabel', group),
                   paste0('longLabel ',group,' tracks for ', nrow(sdf), ' sample x assay combinations'),
                   'group epimap',
                   'subGroup1 view Views ALN=Alignments COV=Coverage',
                   paste(c('subGroup2 sampleType Sample_Type', paste0(smeta$id, '=', smeta$ct)), collapse=' '),
                   paste(c('subGroup3 assayType Assay_Type ', paste0(subassays, '=', subassays)), collapse=' '),
                   'subGroup4 dataType Data_Type Observed=Observed Imputed=Imputed',
                   paste(c('subGroup5 donor Donor ', paste0(dtypes, '=', dtypes)), collapse=' '),
                   'dimensions dimensionX=assayType dimensionY=sampleType',
                   'sortOrder sampleType=+ assayType=+ dataType=+',
                   'dividers sampleType',
                   'dragAndDrop on',
                   'priority 1',
                   'visibility dense',
                   'type bed 3')

    assayview.meta = c(paste0('track CompEpiMapGroup', groupstr, 'ViewCOV'),
                       'shortLabel Coverage',
                       'view COV',
                       'maxHeightPixels 64:32:16',
                       paste0('parent CompEpiMapGroup', groupstr),
                       'visibility dense')

    assay.line = paste(c("\n", assay.meta, "\n"), collapse='\n')
    assayview.line = paste(c("\n", assayview.meta, "\n"), collapse='\n')
    write(assay.line, file=thfile, append=T)
    write(assayview.line, file=thfile, append=T)

    # -----------------------
    # Make subtrack metadata:
    # -----------------------
    # TODO: Do DNase-seq and Histone marks have diff. properties
    # tdf$track = paste0("track BW_",  tdf$id, "_", tdf$mark,"_", substr(tdf$dataset,1,3))
    sdf$track = paste0("track CompGroupBW_",  sdf$id, "_", sdf$mark,"_", substr(sdf$dataset,1,3))
    sdf$subGroups = paste0('subGroups view=COV sampleType=', sdf$id, ' assayType=', sdf$mark, ' dataType=', sdf$dataset, ' donor=', sdf$Donor) 
    sdf$parent = paste0('parent CompEpiMapGroup', groupstr, 'ViewCOV off')
    track.meta = paste(sdf$track, sdf$shortLabel, sdf$longLabel, sdf$bigDataUrl,
                       sdf$attr, sdf$color, sdf$metadata, sdf$subGroups, sdf$parent, sep="\n")

    track.line = paste(track.meta, collapse='\n\n')
    write(track.line, file=thfile, append=T)

}



# -----------------------------
# Make for (new) WUSTL browser:
# -----------------------------
# TODO: Make the metadata definitions:
sdf = tdf
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

thfile = 'wustlhub.fulldataset.config.json'
writeLines(json.line, con=thfile)

grep(",", sdf$ct)
unique(sdf$mark)
unique(sdf$GROUP)
unique(sdf$ct)
unique(sdf$lifestage)
unique(sdf$sex)
unique(sdf$dataset)
unique(sdf$COLOR)
unique(sdf$type)





# -------------------------------
# Additional -  Summary Tracks:
# - Observed - Imputed - CoV?
# -------------------------------




