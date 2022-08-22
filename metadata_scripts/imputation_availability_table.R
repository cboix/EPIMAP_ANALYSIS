#!/usr/bin/R
# Make imputation availability table from ENCODE 
# TODO diff domains
source('~/EPIMAP_ANALYSIS/bin/general_EPIMAP_ANALYSIS.R')
library(reshape)
library(jsonlite)
library(httr)
prefix = paste0(img, '/Imputation/DataSummary_')
rdgy <- colorRampPalette(brewer.pal(n=11,name="RdGy"))(100)

################################################################################################################
consumer_key = "R7GOMSZL";
consumer_secret = "7sccfx4rtdb6sqse";
secret <- RCurl::base64(paste(consumer_key, consumer_secret, sep=":"));
################################################################################################################
# IF DOESNT WORK, GET EXPERIMENTS AS:
# curl -L -u 'R7GOMSZL:7sccfx4rtdb6sqse' -o ~/EPIMAP_ANALYSIS/db/Annotation/all_submitted_released.json https://www.encodeproject.org/search/\?type\=Experiment\&status\=released\&status\=submitted\&replicates.library.biosample.donor.organism.scientific_name\=Homo+sapiens\&format\=json\&limit\=all
################################################################################################################

# ### This is where we obtain most of the metadata
nam = 'all_submitted_released' 
inner = 'type=Experiment&status=released&status=submitted&replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens'
url <- paste("https://www.encodeproject.org/search/?", inner, "&format=json&limit=all", sep="");
# # TODO DOES NOT CURRENTLY WORK to pull data - do in cmdline
# cmd = paste0("curl -L -u 'R7GOMSZL:7sccfx4rtdb6sqse' -o Annotation/", nam, ".json ", url)
# source(cmd)

# # Get and save JSONS for tagAlign and Bam:
for (filetype in c('bam', 'tagAlign', 'bigWig')){
    url_files <- paste0("https://www.encodeproject.org/search/?type=file&status=released&status=submitted&file_format=", filetype, "&frame=object&format=json&limit=all", sep="")
    cmd = paste0("curl -L -u 'R7GOMSZL:7sccfx4rtdb6sqse' -o Annotation/all_files_", filetype, ".json ", url_files)
    print(cmd)
    #system(cmd)
}

# Load in full matrix:
call <- GET(url, config(c("Authorization" = paste("Basic", secret))))
obj <- fromJSON(rawToChar(call$content))
status <- obj$facets[6,"terms"][[1]]

if (status[status$key=='submitted','doc_count'] == 0){ 
    print('Couldn\'t access submitted data, loading previously stored JSON file') 
    obj <- fromJSON(paste0('Annotation/',nam,'.json'))
    status <- obj$facets[6,"terms"][[1]]
}

# Extract important objects:
experiments <- obj[["@graph"]]
experiments$award = experiments$award$project
terms <- obj$facets[[2]]  # Used for extracting the following parts:
# filetypes <- obj$facets[14,"terms"][[1]];

# Keep ChIP-seq ATAC-seq and DNase-seq assays in ENCODE and Roadmap:
assays <- c('ChIP-seq','ATAC-seq', 'DNase-seq', 'total RNA-seq', 'polyA RNA-seq', 'RRBS','WGBS')
experiments <- experiments[experiments$assay_title %in% assays,]
experiments <- experiments[experiments$award %in% c('ENCODE', 'Roadmap'),]

# ====================
# Process experiments:
# ====================
# Marks for prioritization:
sel <- which(experiments$assay_term_name == 'ChIP-seq')
df <- data.frame(id=sel, marks=experiments$target[sel,1])
df$is.mark <- sapply(df$marks, function(x){ length(grep('^H[0-4].*', x)) > 0})
marks <- df$id[df$is.mark]

# Reduce to marks and other:
full <- experiments[c(marks, which(experiments$assay_term_name != 'ChIP-seq')), ]
print(unique(full$target$label))

full <- within(full, assay_target <- assay_title)
full <- within(full, assay_target[assay_title == 'ChIP-seq'] <- target$label[assay_title == 'ChIP-seq'])

# ==============================
# Check which files have bigWig:
# ==============================
f.file='Annotation/all_links_bigWig_ENCODE.rda'
if (length(grep(f.file, dir())) == 1) {
    obj_files <- fromJSON('Annotation/all_files_bigWig.json')
    bwfiles <- obj_files[["@graph"]] 
    bwdata <- bwfiles$dataset
    save(bwfiles, bwdata, file=paste0(f.file))
} else { 
    load(f.file)
}

# Categories: Fold Change, Signal, Both, Any BW?, Nothing

sig = bwfiles$dataset[bwfiles$output_type == 'signal p-value']
fc = bwfiles$dataset[bwfiles$output_type == 'fold change over control']  # subset of signal p-value

full$has.bigWig = full[['@id']] %in% bwfiles$dataset
full$has.signal = full[['@id']] %in% sig
full$bw_file = 'none' 
full$bw_file[full$has.bigWig] = 'other'
full$bw_file[full$has.signal] = 'signal'

# Matrix:
# NOTE:  unclear whether H1-derived and H9-derived neuronal progenitors are distinguished by their biosample term name.
# should we go to biosample summary?

# Add columns to df to generate different overlays: 
# 1) availability
# 2) in Roadmap
# 3) has BigWig for foldchange
# 4) released vs. submitted
# 5) provide approximate storage size
# 6) Priority list, based on: 
# - in ENCODE
# - not just DNase-seq
# - has BigWig files
# - at least 2 (3?) experiments

to.table <- function(df, epord='', ctord='', attribute=''){ 
    formula='accession ~ assay_target + biosample_term_name'
    if (attribute != ''){
        formula = paste(formula, '+', attribute)
    }
    udf <- aggregate(as.formula(formula), length, data=df)
    # Reorder if necessary:
    if (epord == ''){ 
        epdf <- aggregate(accession ~ assay_target, length, data=udf)
        epord <- with(epdf, assay_target[order(accession, decreasing=T)]) 
    }
    if (ctord == ''){ 
        ctdf <- aggregate(accession ~ biosample_term_name, length, data=udf)
        ctord <- with(ctdf, biosample_term_name[order(accession, decreasing=T)]) 
    }
    udf$biosample_term_name <- factor(as.character(udf$biosample_term_name), levels=ctord) 
    udf$assay_target <- factor(as.character(udf$assay_target), levels=epord) 
    udf <- udf[order(udf$biosample_term_name, udf$assay_target),]
    # Reshape to wide:
    if (attribute %in% c('award', 'status', 'bw_file', 'prio')){
        wide <- dcast(udf, assay_target ~ biosample_term_name, value.var=attribute, 
                      fun.aggregate=function(x){ sort(x, decreasing=T)[1]})  # Prioritizes any roadmap expt
    } else if (attribute == '') {
        wide <- dcast(udf, assay_target ~ biosample_term_name, value.var='accession')
    } else {
        wide <- dcast(udf, assay_target ~ biosample_term_name, value.var=attribute)
    }
    mat <- as.matrix(wide[,-1])
    rownames(mat) = wide[,1]
    return(list(mat, epord, ctord))
}

plot.avail <- function(mat, colramp = heat.colors(12), title='', with.rownames=FALSE,
                       add.cutoffs=TRUE, horiz=TRUE, legend=''){
    if (horiz){ 
        mat = t(mat) 
        par(mar=c(4 + with.rownames * 7, 7, 3, 2)) 
    } else {
        par(mar=c(7, 4 + with.rownames * 7, 3, 2)) 
    }
    image(mat, axes=FALSE, main='', col=colramp)
    if (title != '') { mtext(title, side=3, cex=2.25, line=.75)} 
    mtext(paste0('Samples (', ncol(mat) * (!horiz) + nrow(mat) * horiz, ')'), side=2 - horiz * 1, cex=2, line=1.5 + with.rownames * 8)
    mtext(paste0('Genomic Assays (', nrow(mat) * (!horiz) + ncol(mat) * horiz, ')'), side=1 + horiz * 1, cex=2, line=5)
    # grid(nx=nrow(mat), ny=ncol(mat),col='grey',lty='solid',lwd=.25)
    if (horiz){
        text(y=seq(0,1,length.out=ncol(mat)), x=par()$usr[3]+0.01*(par()$usr[4]-par()$usr[3]), labels=colnames(mat), srt=0, adj=1, xpd=TRUE,cex=.75)
    } else { 
        text(x=seq(0,1,length.out=nrow(mat)), y=par()$usr[3]-0.01*(par()$usr[4]-par()$usr[3]), labels=rownames(mat), srt=90, adj=1, xpd=TRUE,cex=.75)
    }
    if (with.rownames * horiz){
        text(x=seq(0,1,length.out=nrow(mat)), y=par()$usr[3]-0.03*(par()$usr[4]-par()$usr[3]), labels=rownames(mat), srt=90, adj=1, xpd=TRUE,cex=.35)
    } else if (with.rownames){
        text(y=seq(0,1,length.out=ncol(mat)), x=par()$usr[3]-0.03*(par()$usr[4]-par()$usr[3]), labels=colnames(mat), srt=0, adj=1, xpd=TRUE,cex=.35)
    }
    if (add.cutoffs){
        if (horiz){ mat <- t(mat) } # revert }
        marks <- rownames(mat)
        N.EPI <- length(marks)
        pct <- rep(0, N.EPI)
        num <- rep(0, N.EPI)
        cuts <- 1:N.EPI
        for (cutoff in cuts){
            mat2 <- mat[marks[1:cutoff],]
            num[cutoff] = sum(!(is.na(mat2)))
            if (cutoff == 1){
                pct[cutoff] = num[cutoff] / length(mat2)
            } else {
                pct[cutoff] = num[cutoff] / prod(dim(mat2))
            }
        }
        breaks = c(7,12,16,19,39)
        bcols <- (breaks - .5)/ (N.EPI-1)
        bpct = pct[breaks]
        bnum = num[breaks]
        ll <- diff(c(0,breaks))
        mids <- ll / 2 + c(0, breaks)[-(length(breaks) +1)]
        textpos <- (mids -.5) / (N.EPI- 1)
        ypos = 0.55 + horiz * 0.2
        projtext = paste0(breaks, ' assays (', round(bpct*100, 1), '%: ', bnum, ')')
        if (horiz) {
            abline(h=bcols,col='darkred',lty='dashed',lwd=2)
            text(y=textpos, x=ypos, labels=projtext, srt=0, adj=0, xpd=TRUE, col=alpha('darkred',.85), cex=2)
        } else { 
            abline(v=bcols,col='darkred',lty='dashed',lwd=2)
            text(x=textpos, y=ypos, labels=projtext, srt=90, adj=0, xpd=TRUE, col=alpha('darkred',.85), cex=2)
        }
    }
    if (legend != ''){
        if (horiz) {
            legend(.1, 1.03, legend=legend, col=colramp, lty=1, lwd=2, cex=1.5, box.lty=0)
        } else {
            legend(.1, 1.03, legend=legend, col=colramp, lty=1, lwd=2, cex=1.5, box.lty=0)
        }
    }
}

core_assays <- c('ChIP-seq','ATAC-seq', 'DNase-seq')
core <- full[full$assay_title %in% core_assays,]
chip <- core[core$assay_title == 'ChIP-seq',]
rdmp <- core[core$award == 'Roadmap',]
enc <- core[core$award == 'ENCODE',]

# Can run for full or just ChIP:
# tag = 'full'
# main = full 
# desc = 'ChIP/DNase/ATAC/BS/RNA-seq'

# tag = 'core'
# main = core
# desc = 'ChIP-seq + DNase/ATAC'

# tag = 'chip'
# main = chip 
# desc  = 'ChIP-seq'

# ===========
# Make plots:
# ===========
tt <- to.table(main)
main.mat <- tt[[1]]

# 1) Availability
pdf(paste0(prefix, tag, '_availability.pdf'), width=20, height=7)
plot.avail(main.mat > 0, col=rdgy[95], with.rownames=TRUE, title=paste(desc, 'Datasets - Availability'))
dev.off()
#plot.avail(main.mat > 0, col=rdgy[95], with.rownames=FALSE)

# 2) In Roadmap
pdf(paste0(prefix, tag, '_award.pdf'), width=20, height=7)
award.mat = to.table(main, epord=tt[[2]], ctord=tt[[3]], attribute='award')[[1]]
award.mat[award.mat == 'ENCODE'] <- 1
award.mat[award.mat == 'Roadmap'] <- 2
class(award.mat) <- 'numeric'
plot.avail(award.mat, col=c('darkblue','firebrick'), with.rownames=TRUE, 
           title=paste(desc,'Datasets - Project'), legend=c('ENCODE','Roadmap'))
dev.off()

# 3) Has BigWig for FoldChange
pdf(paste0(prefix, tag, '_hasBW.pdf'), width=20, height=7)
bw_file.mat = to.table(main, epord=tt[[2]], ctord=tt[[3]], attribute='bw_file')[[1]]
bw_file.mat[bw_file.mat == 'signal'] <- 1
bw_file.mat[bw_file.mat == 'other'] <- 2
bw_file.mat[bw_file.mat == 'none'] <- 3
class(bw_file.mat) <- 'numeric'
plot.avail(bw_file.mat, col=c('darkblue','firebrick', 'orange'), with.rownames=TRUE, add.cutoffs=TRUE, title=paste(desc, 'Datasets - Processed BigWig'),
           legend=c("Signal p-value", "Other (such as read-depth normalized signal)", "No BigWig"))
dev.off()

# 4) Released vs. Submitted
pdf(paste0(prefix, tag, '_status.pdf'), width=20, height=7)
status.mat = to.table(main, epord=tt[[2]], ctord=tt[[3]], attribute='status')[[1]]
status.mat[status.mat == 'submitted'] <- 1
status.mat[status.mat == 'released'] <- 2
class(status.mat) <- 'numeric'
plot.avail(status.mat, col=c('darkblue','firebrick'), with.rownames=TRUE, add.cutoffs=TRUE, 
           title=paste(desc,'Datasets - Release Status'), legend=c('Submitted','Released'))
dev.off()

# 5) Provide approximate storage size
# - Each bw.gz is approximately between 100M and 150M when condensed to 25bp intervals.
# - 1000 experiments (approx roadmap) is between 100G and 150G if stored in bigwig.gz 

# Can run for full or just ChIP:
names <- c('accession', 'assay_title', 'assay_target', 'biosample_term_name', 'award','bw_file', 'status')
df <- full[, names]
tt <- to.table(full)
main.mat <- tt[[1]]
df$assay_target <- factor(as.character(df$assay_target), levels=tt[[2]]) 
df$biosample_term_name <- factor(as.character(df$biosample_term_name), levels=tt[[3]]) 
df <- df[order(df$biosample_term_name, df$assay_target),]

ctnum <- aggregate(accession ~ biosample_term_name, length, data=df)
names(ctnum)[2] <- 'biosample_expt_count'
df <- merge(df, ctnum)

# 6) Priority list, based on: 
# - prioritizing novel ENCODE data followed by Roadmap data
# - celltypes that aren't just DNase-seq
# - experiments that have processed signal tracks
# - celltypes with at least 2 (3?) experiments
df <- within(df, priority <- (biosample_expt_count > 2) + 
             (award == 'ENCODE') +
             (bw_file == 'signal') + 
             (status == 'submitted'))
df <-  df[order(df$priority, decreasing=TRUE),]

# Priority plot:
pdf(paste0(prefix, tag, '_priority.pdf'), width=20, height=7)
df$prio <- paste0('p',df$priority)
priority.mat = to.table(df, epord=tt[[2]], ctord=tt[[3]], attribute='prio')[[1]]
for (i in 1:10){ priority.mat[priority.mat == paste0('p', i)] <- i }
class(priority.mat) <- 'numeric'
num = length(unique(as.numeric(priority.mat))) - 1
lg = sapply(1:num, function(x){paste0(x, ': ', sum(priority.mat == x, na.rm=T), ' files')})
plot.avail(priority.mat, col=c('orange','firebrick','darkblue'), with.rownames=TRUE, title=paste(desc,'Datasets - Priority'), legend=lg)
mtext("(Sum of: 3+ expt. in cell type, is ENCODE dataset, has signal BW file, is submitted)", side=3, cex=1, line=0)
dev.off()

final <- df[,c(names, 'biosample_expt_count', 'priority')]
final <-  final[order(final$priority, decreasing=TRUE),]
write.table(df, file=paste0(prefix, 'ordered_table.tsv'), row.names=F, col.names=F, quote=F, sep="\t")
