#!/usr/bin/R
# Make a full merged DHS set + separate bedfiles for enh/prom

# Load in all relevant data from epimap - enhancers and epigenomes:
# -----------------------------------------------------------------
# Load in all DHS data and merge:
dhsdf = read.delim('masterlist_DHSs_733samples_WM20180608_all_coords_hg19.core.srt.txt', header=F)
srtdf = read.delim('masterlist_DHSs_733samples_WM20180608_all_coords_hg19_r25_e100_names.core.srt.tsv', header=T)
mapdf = read.delim('masterlist_DHSs_733samples_WM20180608_all_chunkIDs2indexIDs.txt', header=F)
names(dhsdf) = c('chr','start','end','name')
names(mapdf) = c('id','name')

df = merge(srtdf, dhsdf, all.x=TRUE)
df = merge(df , mapdf, all.x=TRUE)
df = df[order(df$cls), c('chr','start','end','cls','name','id')]
names(df)[4] = 'row'
df$strand = '.'

# Load in enhancer locations + other sets and add them to the full dataframe:
# ---------------------------------------------------------------------------
eind = as.numeric(scan('ENH_masterlist_indices_0indexed.tsv','c')) + 1
pind = as.numeric(scan('PROM_masterlist_indices_0indexed.tsv','c')) + 1
dind = as.numeric(scan('DYADIC_masterlist_indices_0indexed.tsv','c')) + 1

df$is.enh = 0
df$is.prom = 0
df$is.dyadic = 0
df$is.enh[eind] = 1
df$is.prom[pind] = 1
df$is.dyadic[dind] = 1

# Save outputs:
# -------------
# Write this master table out:
write.table(df, 'masterlist_DHSs_733samples_WM20180608_all_coords_hg19_epimap_annotated.tsv', quote=F, row.names=F, col.names=T, sep="\t")

# Make BED files for each element type:
cols = c('chr','start','end','strand', 'name', 'row')
edf = df[df$is.enh == 1, cols]
pdf = df[df$is.prom == 1, cols]
ddf = df[df$is.dyadic == 1, cols]

write.table(edf, 'ENH_masterlist_locations.bed', quote=F, row.names=F, col.names=F, sep="\t")
write.table(pdf, 'PROM_masterlist_locations.bed', quote=F, row.names=F, col.names=F, sep="\t")
write.table(ddf, 'DYADIC_masterlist_locations.bed', quote=F, row.names=F, col.names=F, sep="\t")

