#!/usr/bin/R 
# Server side for imputation visualization:
library(shiny)
options(scipen=45)
domain <- Sys.getenv('DOMAINNAME')
if (domain == 'broadinstitute.org'){ 
    img <- '/broad/compbio/cboix/EPIMAP_ANALYSIS/img/'
    bin <- '/broad/compbio/cboix/EPIMAP_ANALYSIS/bin/'
} else { 
    img <- '~/EPIMAP_ANALYSIS/img/'
    bin <- '~/EPIMAP_ANALYSIS/bin/'
}

# Load functions:
source('imputation_avail_funcs.R')

# Load distance matrix
mat = read.delim(paste0('data/DataSummary_ordered_distances.tsv'), header=T)
mat = as.matrix(mat)
colnames(mat) = as.character(unlist(read.delim(paste0('data/celltypes.tsv'), header=F)))
rownames(mat) = colnames(mat)
dt = dist(1 - mat, 'cosine') 
ht <- hclust(dt, 'ward.D2')
cocl <- order.optimal(dt, ht$merge)$order 

# Load annotations:
df = read.delim(paste0('data/DataSummary_ordered_table.tsv'), header=F, stringsAsFactors=F)
names(df) <- c('biosample_term_name','accession','assay_title','assay_target','award', 'bw_file','status','biosample_expt_count','priority','prio')
ct = as.character(unique(df$biosample_term_name))
core_assays <- c('ChIP-seq','ATAC-seq', 'DNase-seq', 'polyA RNA-seq', 'total RNA-seq')
df <- df[df$assay_title %in% core_assays,]

# Keep direct matches between two tissues:
newanno <- read.delim('data/new_samples_mapping.tsv',header=T)
newanno <- newanno[as.character(newanno$Snyder.samples) != '',]
newanno <- newanno[as.character(newanno$ENTEx) != '',]

# Load in namesets
core = c('IMR-90', 'H1-hESC','H9', 'trophoblast cell', 'mesenchymal stem cell', 'neural stem progenitor cell', 'K562', 'adrenal gland', 'heart left ventricle')
inc = c('upper lobe of left lung','cerebellum','HEK293','kidney epithelial cell','mononuclear cell','Caco-2','cardiac mesoderm','ES-I3','GM06990','BE2C','prostate','vagina','neuroepithelial stem cell')
newsamp = lapply(newanno$Snyder, function(x){which(ct %in% x)})
matches = unlist(newsamp)
newct = ct[matches]
divlist = c()
# potential = c('BE2C','Caco-2','ES-I3','HL-60','LNCaP clone FGC','KMS-11','NCI-H929','SK-N-MC','GM06990','WERI-Rb-1','RWPE2','MG63','SJCRH30','SJSA1','HAP-1','DND-41','OCI-LY7','uterus','vagina','ovary','omental fat pad','astrocyte','brain microvascular endothelial cell','cardiac fibroblast','fibroblast of lung')
potential = c('BE2C','Caco-2','ES-I3','HL-60','KMS-11','NCI-H929','GM06990','WERI-Rb-1','RWPE2','MG63','SJCRH30','SJSA1','HAP-1','DND-41','OCI-LY7','uterus','vagina','omental fat pad','brain microvascular endothelial cell','cardiac fibroblast')


shinyServer(function(input, output,session) {
                values <- reactive({
                    cells <- sort(unique(df$biosample_term_name))
                    cols = c('biosample_term_name','assay_target')
                    subdf = df[df$biosample_term_name %in% cells,]
                    size.table = data.frame(step=0, ndata=nrow(subdf), ncomb=nrow(unique(subdf[,cols])),
                                            ncell=length(cells), nassay=length(unique(subdf$assay_target)))
                    # Starting from full dataset (df, mat), prune: 
                    # 1. At least 2 assays
                    if (input$gt2) {
                        cdf = aggregate(assay_target ~ biosample_term_name, subdf,
                                        function(x){length(unique(x))})
                        cells <- cdf$biosample_term_name[cdf$assay_target > 1] 
                        subdf = df[df$biosample_term_name %in% cells,]
                        size.table <- rbind(size.table, c(1, nrow(subdf), nrow(unique(subdf[,cols])), 
                                                          length(cells), length(unique(subdf$assay_target))))
                    }
                    # 2. Remove datasets not missing core targets
                    if (input$miss){
                        marks=c('DNase-seq','H3K4me1','H3K4me3','H3K36me3','H3K27me3','H3K9me3','H3K27ac')
                        mdf = subdf[subdf$assay_target %in% marks,]
                        cdf = aggregate(assay_target ~ biosample_term_name, mdf,
                                        function(x){length(unique(x))})
                        # Keep cells with at least one core mark:
                        cells = cdf$biosample_term_name[cdf$assay_target < length(marks)]
                        subdf = df[df$biosample_term_name %in% cells,]
                        size.table <- rbind(size.table, c(2, nrow(subdf), nrow(unique(subdf[,cols])), 
                                                          length(cells), length(unique(subdf$assay_target))))
                    }
                    # 3. Keep only those with signal tracks
                    if (input$sign){
                        rows = which(!(subdf$assay_title == 'ChIP-seq' &
                                       subdf$bw_file != 'signal'))
                        subdf = subdf[rows,]
                        cells = sort(unique(subdf$biosample_term_name))
                        size.table <- rbind(size.table, c(3, nrow(subdf), nrow(unique(subdf[,cols])), 
                                                          length(cells), length(unique(subdf$assay_target))))
                    }
                    set <- input$txt
                    finaldf = df[df$biosample_term_name %in% set,]
                    size.table <- rbind(size.table, c(4, nrow(finaldf), nrow(unique(finaldf[,cols])),
                                                      length(set), length(unique(finaldf$assay_target))))
                    if (input$core){ set <- c(set, core) }
                    finaldf = df[df$biosample_term_name %in% set,]
                    size.table <- rbind(size.table, c(5, nrow(finaldf), nrow(unique(finaldf[,cols])), 
                                                      length(set), length(unique(finaldf$assay_target))))
                    if (input$sugg){ set <- c(set, potential) }
                    finaldf = df[df$biosample_term_name %in% set,]
                    size.table <- rbind(size.table, c(6, nrow(finaldf), nrow(unique(finaldf[,cols])), 
                                                      length(set), length(unique(finaldf$assay_target))))
                    set = unique(set)
                    rn = rownames(mat)
                    cells = rn[rn %in% c(set, cells)]
                    # Pick up to N diverse.
                    set = pickToN(mat[cells,cells], set, input$diverse)
                    set = unique(set)
                    print(core[core %in% set])
                    finaldf = df[df$biosample_term_name %in% set,]
                    size.table <- rbind(size.table, c(7, nrow(finaldf), nrow(unique(finaldf[,cols])), 
                                                      length(set), length(unique(finaldf$assay_target))))
                    print(size.table)
                    return(set)
                })

                # Get matrix for availability plots:
                processList <- function(values, inputAttr, inputSignal){
                    coredf = df[df$biosample_term_name %in% values,]
                    tt <- to.table(coredf)
                    if (inputSignal){
                        rows = which(!(coredf$assay_title == 'ChIP-seq' &
                                       coredf$bw_file != 'signal'))
                        coredf = coredf[rows,]
                    }
                    wide.mat = tt[[1]]
                    wide.mat = wide.mat > 0
                    if (inputAttr != 'accession'){
                        wide.mat = to.table(coredf, epord=tt[[2]], ctord=tt[[3]],
                                            attribute=inputAttr)[[1]]
                    }
                    return(list(wide.mat, coredf))
                }
                matlist <- reactive({
                    return(processList(values(), input$attr, input$sign))
                })

                # Get matrix for availability plots:
                getTable <- reactive({ 
                    coredf = df[df$biosample_term_name %in% values(),]
                    if (input$sign){
                        rows = which(!(coredf$assay_title == 'ChIP-seq' &
                                       coredf$bw_file != 'signal'))
                        coredf = coredf[rows,] }
                    return(coredf)
                })

                # Get matrix for availability plots:
                getBreaks <- reactive({
                    return(get.breaks(ht, cocl, input$nbreaks))
                })

                processAvailMat <- function(inputTitle, inputAttr, wide.mat){
                    if (inputAttr == 'accession'){
                        colors = rdgy[95]
                        legend = ""
                        title = paste(inputTitle, 'Datasets - Availability')
                    } else if (inputAttr == 'status') {
                        wide.mat[wide.mat == 'submitted'] <- 1
                        wide.mat[wide.mat == 'released'] <- 2
                        colors = c('darkblue','firebrick')
                        legend = c('Submitted','Released')
                        title = paste(inputTitle, 'Datasets - Release Status')
                    } else if (inputAttr == 'award') {
                        wide.mat[wide.mat == 'ENCODE'] <- 1
                        wide.mat[wide.mat == 'Roadmap'] <- 2
                        colors = c('darkblue','firebrick')
                        legend = c('ENCODE','Roadmap')
                        title = paste(inputTitle, 'Datasets - Project')
                    } else if (inputAttr == 'bw_file') {
                        wide.mat[wide.mat == 'signal'] <- 1
                        wide.mat[wide.mat == 'other'] <- 2
                        wide.mat[wide.mat == 'none'] <- 3
                        colors = c('darkblue','firebrick', 'orange')
                        legend = c("Signal p-value",
                                   "Other (such as read-depth normalized signal)", "No BigWig")
                        title = paste(inputTitle, 'Datasets - Processed BigWig')
                    } else if (inputAttr == 'priority') { } # TODO Priority figure
                    class(wide.mat) <- 'numeric'
                    return(list(wide.mat, colors, legend, title))
                }
                availAttr <- reactive({
                    return(processAvailMat(input$title, input$attr, matlist()[[1]]))
                })

                #==========================================
                #-------------- PLOTS ---------------------
                #==========================================
                output$Avail <- renderPlot({
                    vals <- values()
                    attrlist <- availAttr()
                    wide.mat = attrlist[[1]]
                    colors = attrlist[[2]]
                    legend = attrlist[[3]]
                    title = attrlist[[4]]
                    plot.avail(wide.mat, col=colors, with.rownames=TRUE,
                               title=title, legend=legend, cex.lab='auto', 
                               highlightct=potential, 
                               highlightct2=core)
                },height=800,width=1500)


                # TODO is it possible to only refresh the axis labels?
                output$Heatmap <- renderPlot({
                    vals <- values()
                    breaks <- getBreaks()
                    plot.sym(mat, colramp=spec, breaks=breaks, yaxlab=colnames(mat),
                             quant=0.02, highlight=vals)
                    title(main=paste(input$title, 'datasets within distance matrix'))
                },height=1400,width=1700)

                output$Table = renderDataTable(getTable())

                output$downAvail <- downloadHandler(filename = 'availability.pdf',
                                                    content = function(file){
                                                        pdf(file, width=15, height=8)
                                                        vals <- values()
                                                        attrlist <- availAttr()
                                                        wide.mat = attrlist[[1]]
                                                        colors = attrlist[[2]]
                                                        legend = attrlist[[3]]
                                                        title = attrlist[[4]]
                                                        plot.avail(wide.mat, col=colors, with.rownames=TRUE,
                                                                   title=title, legend=legend, cex.lab='auto',
                                                                   highlightct=potential, highlightct2=core)
                                                        dev.off()
                                                    })

                output$downHeat <- downloadHandler(filename = 'heatmap.pdf',
                                                    content = function(file){
                                                        pdf(file, width=18, height=15)
                                                        vals <- values()
                                                        breaks <- getBreaks()
                                                        plot.sym(mat, colramp=spec, breaks=breaks, yaxlab=colnames(mat),
                                                                 quant=0.02, highlight=vals)
                                                        title(main=paste(input$title, 'datasets within distance matrix'))
                                                        dev.off()
                                                    })

                output$downTable <- downloadHandler(filename = 'availability_table.tsv',
                                                    content = function(file){
                                                        write.table(getTable(), file, row.names=F,
                                                                    col.names=T, quote=F, sep="\t")
                                                    })

})
