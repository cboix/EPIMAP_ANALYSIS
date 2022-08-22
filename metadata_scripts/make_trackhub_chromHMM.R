#!/usr/bin/R
# ---------------------------
# Make the ChromHMM trackhub:
# ---------------------------
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


# -------------------------
# Reduce to 833 epigenomes:
# -------------------------
agemap = rd.map[,c('BSSID','age','ageUnits')]
names(agemap)[1] = 'id'
meta = merge(meta, agemap, all.x=TRUE)
rownames(meta) = meta$id
meta = meta[sort(cellorder),]
meta$fullage = paste0(meta$age, ' ', meta$ageUnits, 's')
meta$fullage[grep('unknown',meta$fullage)] = 'unknown'
meta$ctspace = gsub("_", " ", meta$ct)
meta$ctspace = gsub("  ", " ", meta$ctspace)

for (genome in c('hg19','hg38')){
    # Making single file:
    hub.meta = c('hub ChromHMM',
                 'shortLabel EpiMap release v0.9 (ChromHMM)',
                 'longLabel ChromHMM tracks for EpiMap v0.9',
                 'useOneFile on',
                 'email cboix@mit.edu',
                 '\n', 
                 paste('genome', genome))

    chmmdir = paste0('https://personal.broadinstitute.org/cboix/epimap/ChromHMM/observed_aux_18_', genome, '/CALLS/bigBed/')

    # Make subtrack metadata:
    tdf = data.frame(meta$id)
    tdf$track = paste0("track ChromHMM_", meta$id)
    tdf$shortLabel = paste("shortLabel", meta$infoline)
    # tdf$longLabel = paste("longLabel", meta$id, meta$ctspace)
    tdf$longLabel = paste("longLabel", meta$id, meta$name)  # TODO: test if works
    tdf$bigDataUrl = paste0('bigDataUrl ', chmmdir, meta$id, '_18_CALLS_segments.bb')
    tdf$attr = paste(c("type bigBed 9 .", 
                       "itemRgb on",
                       "visibility dense", 
                       "noInherit on", 
                       "maxHeightPixels 64:32:16"), collapse='\n')
    tdf$metadata = paste('metadata', 
                         paste0('"Epigenome ShortName"="', meta$infoline,'"'),
                         paste0('"Epigenome FullName"="', meta$name,'"'),
                         paste0('"Group"="<span style="color:', meta$COLOR,'">', meta$GROUP, '</span>"'),
                         paste0('"Type"="', meta$type,'"'),
                         paste0('"Sex"="', meta$sex,'"'),
                         paste0('"Age"="', meta$fullage,'"'),
                         paste0('"Donor"="', meta$Donor,'"'),
                         sep=" ")
    tdf$group = as.character(meta$GROUP)
    tdf$id = as.character(meta$id)
    tdf$ct = as.character(meta$ct)
    tdf$Donor = as.character(meta$Donor)
    tdf$subGroups = paste0('subGroups view=ChromHMM sampleType=', meta$id) 
    tdf$parent = paste0('parent EpiMapConsolidatedHMMViewHMM off')
    track.meta = paste(tdf$track, tdf$shortLabel, tdf$longLabel, tdf$bigDataUrl,
                       tdf$attr, tdf$metadata, tdf$subGroups, tdf$parent, sep="\n")

    # Overall metadata:
    ovl.meta = c('track EpiMapConsolidatedHMM',
          'compositeTrack on',
          'shortLabel chromHMM',
          'longLabel chromHMM tracks from EpiMap',
          'subGroup1 view Views ChromHMM=ChromHMM',
          paste(c('subGroup2 sampleType Sample_Type', paste0(meta$id, '=', meta$ct)), collapse=' '),
          'dimensions dimensionX=view dimensionY=sampleType',
          'sortOrder sampleType=+ view=+',
          'dividers sampleType',
          'dragAndDrop on',
          'visibility hide',
          'type bed 3')

    # View track:
    view.meta = c('track EpiMapConsolidatedHMMViewHMM',
                  'shortLabel ChromHMM',
                  'view ChromHMM',
                  'maxHeightPixels 64:32:16',
                  'parent EpiMapConsolidatedHMM',
                  'visibility dense')

    # --------------------------------
    # Write out the chromHMM trackHub:
    # --------------------------------
    thfile = paste0('trackHub_chromHMM_', genome, '.txt')

    hub.line = paste(hub.meta, collapse='\n')
    ovl.line = paste(ovl.meta, collapse='\n')
    view.line = paste(view.meta, collapse='\n')
    track.line = paste(track.meta, collapse='\n\n')

    writeLines(c(hub.line, ovl.line, view.line, track.line), con=thfile, sep = "\n\n")

    # ---------------------
    # Additional -  ByGroup
    # ---------------------
    # By Assay:
    super.meta = c('track CompEpiMapbyGroup',
                   'superTrack on',
                   'shortLabel By Group',
                   'longLabel EpiMap ChromHMM by Group',
                   'group epimap',
                   'visibility dense',
                   'priority 11')

    super.line = paste(c("\n", super.meta, "\n"), collapse='\n')
    write(super.line, file=thfile, append=T)

    groups = sort(unique(as.character(meta$GROUP)))
    for (group in groups){
        print(group)
        astr = sub(" & ", "and", group)
        astr = sub("\\. ", "", astr)
        sdf = tdf[tdf$group == group,]
        smeta = unique(sdf[,c('id','ct')])

        ddf = aggregate(id ~ Donor, sdf, length)
        ddf = ddf[order(ddf$id, decreasing=T),]
        dtypes = ddf$Donor

        # Overall metadata:
        group.meta = c(paste0('track CompEpiMap', astr),
                       'compositeTrack on',
                       'parent CompEpiMapbyGroup',
                       paste('shortLabel', group),
                       paste('longLabel',group,'tracks for', nrow(smeta), 'sample type(s)'),
                       'subGroup1 view Views ChromHMM=ChromHMM',
                       paste(c('subGroup2 sampleType Sample_Type', paste0(smeta$id, '=', smeta$ct)), collapse=' '),
                       'dimensions dimensionX=view dimensionY=sampleType',
                       'sortOrder sampleType=+ view=+',
                       'dividers sampleType',
                       'priority 1',
                       'dragAndDrop on',
                       # 'type bigWig 0 1.0',
                       'type bed 3',
                       # "type bigBed 9 .", 
                       'visibility dense')

        # View track:
        view.meta = c(paste0('track CompEpiMap', astr, 'ViewHMM'),
                      paste('shortLabel ChromHMM for', group),
                      'view ChromHMM',
                      'maxHeightPixels 64:32:16',
                      paste0('parent CompEpiMap', astr),
                      'visibility dense')

        group.line = paste(c("\n", group.meta, "\n"), collapse='\n')
        view.line = paste(c("\n", view.meta, "\n"), collapse='\n')
        write(group.line, file=thfile, append=T)
        write(view.line, file=thfile, append=T)

        # Make subtrack metadata:
        sdf$track = paste0("track CompGroupHMM_", astr, '_',  sdf$id)
        sdf$parent = paste0('parent CompEpiMap', astr, 'ViewHMM on')
        track.meta = paste(sdf$track, sdf$shortLabel, sdf$longLabel, sdf$bigDataUrl,
                           sdf$attr, sdf$metadata, sdf$subGroups, sdf$parent, sep="\n")

        track.line = paste(track.meta, collapse='\n\n')
        write(track.line, file=thfile, append=T)
    }
}

