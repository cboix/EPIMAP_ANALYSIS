#!/usr/bin/R
source('~/EPIMAP_ANALYSIS/bin/general_EPIMAP_ANALYSIS.R')

meta <- read.delim('expression/metadata.tsv',header=T,sep="\t")
meta <- meta[meta$File.Status == 'released',]

meta$UID <- paste(meta$Experiment.accession,
                  meta$Biological.replicate.s.,
                  meta$Technical.replicate,
                  sep=".")

# Fields we are interested in:
nam <- c('UID',
         'File.accession',
         'File.format',
         'Output.type',
         'Experiment.accession',
         'Biosample.term.id',
         'Biosample.term.name',
         'Biosample.subcellular.fraction.term.name',
         'Library.made.from', # Important for RNA, unless only using gene quantifications.
         'Library.depleted.in',
         'Derived.from',
         "Biological.replicate.s.",
         "Technical.replicate",
         'Assembly',
         'File.Status',
         'File.download.URL',
         'Lab')
dat <- meta[,nam]

# File types:
ft.df <- unique(dat[,c('File.format','Output.type')])
print(ft.df)

ft.of.interest <- c('gtf','tsv')
dat <- dat[dat$File.format %in% ft.of.interest,]

bioM <- unique(meta[,c('Biosample.term.id','Biosample.term.name')])
bio.df <- unique(dat[,c('Biosample.term.id','Biosample.term.name')])

# ONLY GENE QUANTIFICATIONS:
tx <- meta[meta$Output.type == 'gene quantifications',nam]
bioTX <- unique(tx[,c('Biosample.term.id','Biosample.term.name')]) # Have 179 quantifications...

# Uniquely identified from EXPERIMENT:
uids <- unique(tx$UID)

tx$Assembly <- factor(tx$Assembly,levels=c('GRCh38','hg19'))
tx <- tx[order(tx$Assembly),]

eids <- c()
for (uid in uids){ 
    ids <- which(tx$UID == uid)
    if (length(ids) == 1){
        eids <- c(eids,ids)
    } else {
        if (length(ids) == 4){ 
            eids <- c(eids,ids[1:2])
        } else { 
            eids <- c(eids,ids[1])
        }
    }
}

exp.df <- tx[eids,]

# SAVE REFS:
write(as.character(exp.df$File.download.URL),'expression/TX_links.tsv')

# SAVE FINAL INFO TABLE:
write.table(exp.df,'expression/TX_all_submitted_released.tsv',row.names=F,quote=F,sep="\t")

