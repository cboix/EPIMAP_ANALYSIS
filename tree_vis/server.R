#!/usr/bin/R 
# -----------------------------------------
# Server side for GWAS - tree visualization
# -----------------------------------------
library(shiny)
options(scipen=45)
domain <- Sys.getenv('DOMAINNAME')
if (domain == 'broadinstitute.org'){
    maindir = '/broad/compbio/cboix/EPIMAP_ANALYSIS/'
} else {
    maindir = '~/EPIMAP_ANALYSIS/'
}
dbdir  <- paste0(maindir, 'db/')
img <- paste0(maindir, 'img/')
bin <- paste0(maindir, 'bin/')
pdfdir = paste0(img, "gwas_tree_analysis/examples/")
# pdfpref = paste0(pdfdir, 'enhancers_e2500_cons_parent_')
addResourcePath('pdf_folder', pdfdir)
pdfpref = paste0('pdf_folder/', 'enhancers_e2500_cons_parent_')

# --------------------------------------------------------
# Load in + process all of the relevant matrices/datasets:
# --------------------------------------------------------
usetree = 'enhancers'
tol = 2500
singlematch = FALSE
plotting.only = TRUE
args = c(usetree, tol, singlematch, plotting.only)
commandArgs <- function(trailingOnly=TRUE){
    return(c(usetree, tol, singlematch, plotting.only)) }
print(commandArgs())
source(paste0(bin, 'load_gwastree_analysis.R'))

# ----------------------------------
# Load the pre-computed enrichments:
# ----------------------------------
type = 'cons'
against = 'parent'
apref = paste0(type, '_', against)
etfile = paste0(dbdir, 'logreg_', apref, midpref, 'pvals_elist.Rda')
load(etfile)

snpuidlist = scan(paste0(dbdir,'full_snpuidlist.txt'), 'c', sep="\n")

# Other functions:
setup.circleplot = function(size=1){
    # Set up plotting region:
    plot.new()
    circle_size = unit(size, "snpc") # snpc unit gives you a square region
    pushViewport(viewport(x=0, y=0.5, width=circle_size, height=circle_size, just=c("left", "center")))
    par(omi = gridOMI(), new = TRUE)
    return(circle_size)
}

shinyServer(function(input, output,session) {
                # Dynamically get trait data:
                trait.ll <- reactive({
                    print(input$trait)
                    return(lrlist[[input$trait]])
                })

                # Dynamically update dendrogram:
                dend.dd <- reactive({
                    return(setup_dendrogram(dend3, trait.ll(), udf, declist=declist, 
                                            bcutoff=3, altline=(3 * (!input$branch))))
                })

                # Dynamically process other options for plotting:
                getopts <- reactive({
                    ll = trait.ll()
                    if (input$leafpval){ leaf.pval = leafmat[input$trait,] } else { leaf.pval = NULL }
                    title = paste(input$trait, '-', ll$title)
                    opts = list(leaf.pval=leaf.pval,
                                fulltrack=input$fulltrack,
                                meta=input$meta,
                                legend=input$legend,
                                addtree=input$addtree, 
                                title=title)
                    return(opts)
                })

                get.pdfopts <- reactive({
                    ind = which(snpuidlist == input$uid)
                    ipref = sprintf("%05d", ind)
                    # pdftail = paste0(ipref, '_lrtree')
                    pdftail = paste0(ipref, '_lrtree_adj')
                    if (input$allgenes){
                        pdftail = paste0(pdftail,'_genes_all')
                    } else {
                        if (input$genes){ pdftail = paste0(pdftail, '_genes') }
                        if (input$minp){ pdftail = paste0(pdftail, '2') }
                    }
                    pdftail = paste0(pdftail, '.pdf')
                    return(paste0(pdfpref, pdftail))
                })

                #==========================================
                #-------------- PLOTS ---------------------
                #==========================================
                # Most of the computation happens in the plotting.
                # TODO: Find a way to speed up/cache the plotting
                output$Tree <- renderPlot({
                    opts <- getopts()
                    dd = dend.dd()
                    # Make plot, with options:
                    circle_size = setup.circleplot()
                    circleplot(dend=dd$dend, lab=lab, udf=dd$udf, 
                               fractional=FALSE, fulltrack=opts$fulltrack,
                               add.metadata=opts$meta, with.tree=opts$addtree,
                               hit.track=dd$hitdf, leaf.pval=opts$leaf.pval)
                    upViewport()
                    circos.clear()
                    if (opts$legend){
                        draw(pd.legend.ext, x = circle_size, just = "left")
                    }
                    title(opts$title)
                }, height=1200, width=1450)

                output$downTree <- downloadHandler(filename = 'gwastree.pdf',
                                                   content = function(file){
                                                       pdf(file, width=14.5, height=12)
                                                       opts <- getopts()
                                                       dd = dend.dd()
                                                       # Make plot, with options:
                                                       circle_size = setup.circleplot()
                                                       circleplot(dend=dd$dend, lab=lab, udf=dd$udf, 
                                                                  fractional=FALSE, fulltrack=opts$fulltrack,
                                                                  add.metadata=opts$meta, with.tree=opts$addtree,
                                                                  hit.track=dd$hitdf, leaf.pval=opts$leaf.pval)
                                                       upViewport()
                                                       circos.clear()
                                                       if (opts$legend){
                                                           draw(pd.legend.ext, x = circle_size, just = "left")
                                                       }
                                                       title(opts$title)
                                                       dev.off()
                                                   })


                output$Pdf <- renderUI({
                    pdf.file <- get.pdfopts()
                    print(pdf.file)
                    # tags$iframe(src=pdf.file)
                    tags$iframe(src=pdf.file, width="1450", height="1200")
                })
                # }, height=1200, width=1450)

})
