#!/usr/bin/R
# UI side for imputation visualization:
library(shiny)
library(shinyjs)
library(shinyTree)
library(shinythemes)
library(shinydashboard)
library(shinycssloaders)
library(shinyWidgets)
library(visNetwork)
library(markdown)
library(DT)
library(heatmaply)
options(scipen=45)

# ----------------
# For sample tree:
# ----------------
# Make types list:
load('data/metadata_app.Rda')
meta$infostr = paste0(meta$id, ': ', meta$info)
meta$ctstr = paste0(meta$id, ': ', meta$ct)
gcols = colvals$group
groups.full = names(gcols)
groups = gsub("[& -.]","", tolower(groups.full))
names(gcols) = groups
typeslist = c("{", sapply(groups, function(x){ paste0("'",x,"': {'a_attr' : { 'style' : 'color:", gcols[x],"' }}, ") }), "}")
typeslist = paste(typeslist, collapse="")
co859 = scan('data/859_samples.txt',quiet=TRUE, 'c', sep="\n")
cellorder = scan('data/833_samples.txt',quiet=TRUE, 'c', sep="\n")
meta = meta[meta$id %in% co859,]

# For picking marks:
all.assays = scan('data/all_assays.txt', quiet=TRUE, 'c', sep="\n")

# -------------------------
# For the preset trackhubs:
# -------------------------
resdir = "http://personal.broadinstitute.org/cboix/epimap/"
UCSCprefix = 'http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&hubUrl=http://personal.broadinstitute.org/cboix/epimap/trackhubs/'
UCSC38prefix = 'http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&hubUrl=http://personal.broadinstitute.org/cboix/epimap/trackhubs/'
WUSTLprefix = 'http://epigenomegateway.wustl.edu/legacy/?genome=hg19&datahub=http://personal.broadinstitute.org/cboix/epimap/trackhubs/'
WUSTLnewprefix = 'http://epigenomegateway.wustl.edu/browser/?genome=hg19&hub=http://personal.broadinstitute.org/cboix/epimap/trackhubs/'

# -------------------------
# For the GWAS enrichments:
# -------------------------
# Read in autocomplete list of GWAS:
# Read in the list of annotations to look at:
tree.gwas = scan('data/tree_gwas_list.txt', quiet=TRUE, 'c', sep="\n")
tree.gwas1p = scan('data/tree_gwas1p_list.txt', quiet=TRUE, 'c', sep="\n")

# For checkbox input in grid-style
tweaks <- list(tags$head(tags$style(HTML("
                                         .multicol { 
                                             height: 350px;
                                             -webkit-column-count: 3; /* Chrome, Safari, Opera */ 
                                                 -moz-column-count: 3;    /* Firefox */ 
                                                 column-count: 3; 
                                             -moz-column-fill: auto;
                                             -column-fill: auto;
                                         } 
                                         ")) 
                                         ))
group.checkboxes <- tags$div(align = 'left', class = 'multicol', 
                             checkboxGroupInput(inputId  = 'selGroups', 
                                                label    = "Restrict to sample groups:", 
                                                choices  = groups.full,
                                                selected = NULL,
                                                inline   = FALSE))

# UI for tree editor:
shinyUI(
        # fluidPage(
        navbarPage(theme = shinytheme("cosmo"),
                   "Epimap Data Exploration", id='navBarSet',
                   tabPanel("Epigenomes Table + TrackHubs", value='tracksNV',
                            sidebarLayout(sidebarPanel(
                                                       selectizeInput('textNodes', 'Add specific samples (by short name):', choices = meta$infostr, multiple = TRUE, options = list(create = TRUE)),
                                                       selectizeInput('textNodesCT', 'Add specific samples (by full name):', choices = meta$ctstr, multiple = TRUE, options = list(create = TRUE)),
                                                       radioButtons("selTree", "Select from:", c("Groups" = "groups", "Tree" = "tree"), inline=TRUE, selected='groups'),
                                                       tags$b("Add samples from sample groups or", tags$a(href=paste0(resdir, 'gwas_resources/enhancer_sample_tree.png'), "enhancer tree:")),
                                                       shinyTree("tree", "Add samples from tree:", checkbox=TRUE, theme="proton", types=typeslist)),
                                          mainPanel(tabsetPanel(type='tabs',
                                                                tabPanel("Sample Table", 
                                                                         DT::dataTableOutput("selTable"),
                                                                         downloadButton('downSampleTable', label="Download Table")),
                                                                tabPanel("Preset TrackHubs and Views", 
                                                                         h4("Representative track hubs:"),
                                                                         tags$b('1. Full track hubs (all 17,938 tracks):'),
                                                                         # Full track hub here:
                                                                         tags$a(href=paste0(WUSTLnewprefix, 'wustlhub.fulldataset.config.json'), "WUSTL"),'/',
                                                                         tags$a(href=paste0(UCSCprefix, 'trackHub_full.txt'), "UCSC"),
                                                                         br(),
                                                                         "These full hubs are huge track hubs and will be very slow to load.", br(),
                                                                         "We recommend making custom hubs (next tab) to restrict to datasets of interest.",
                                                                         br(), br(),

                                                                         tags$b('2. Average tracks for 33 major sample groups in each mark:'),
                                                                         tags$a(href=paste0(WUSTLprefix, 'wustlhub.groupaverages.json'), "WUSTL"),'/',
                                                                         tags$a(href=paste0(UCSCprefix, 'trackHub_groupaverages.txt'), "UCSC"),
                                                                         br(), br(),
                                                                         tags$b('3. Chromatin state tracks:'), tags$a(href="https://epilogos.altius.org/?application=viewer&sampleSet=vC&mode=single&genome=hg19&model=18&complexity=KL&group=all&chrLeft=chr5&chrRight=chr5&start=38580123&stop=39602208", "Epilogos"),'/',
                                                                         # tags$a(href=paste0(WUSTLprefix, 'wustlhub.chromhmm.json'), "WUSTL"),'/',
                                                                         tags$a(href=paste0(UCSCprefix, 'trackHub_chromHMM_hg19.txt'), "UCSC (hg19)"),'/',
                                                                         tags$a(href=paste0(UCSC38prefix, 'trackHub_chromHMM_hg38.txt'), "UCSC (hg38)"),
                                                                         br(), br(),
                                                                         tags$b('4. All epigenomes and average tracks for a specific sample group:'), 
                                                                         uiOutput("presetHubGroup"),
                                                                         selectizeInput('selTrackHubGroups', 'Selected sample group:',
                                                                                        choices = groups.full, multiple = FALSE, options = list(create = TRUE))

                                                                ),
                                                                tabPanel("Make Custom TrackHubs", 
                                                                         # Choose marks here:
                                                                         h3("TrackHub Options:"),
                                                                         # Only imputed, only observed, both
                                                                         fluidRow(
                                                                         column(2,radioButtons("selIO", "Track type(s):",
                                                                                      c("Imputed" = "imp", "Observed" = "obs", "Both" = "both"), selected='both')),
                                                                         column(2, checkboxGroupInput("selMarkTiers", "Marks/Assays (by set):",
                                                                                            c("Tier 1 (Core Marks + DNase-seq)" = "t1", "Tier 2 (Secondary Marks + ATAC-seq)" = "t2", "Tier 3 (DNA Factors)" = "t3", "Tier 4 (Other Marks)" = "t4"),selected='t1')),
                                                                         column(4, selectizeInput('selMarks', 'Specific histone marks and assays to add to the track hub:',
                                                                                        choices = all.assays, multiple = TRUE, options = list(create = TRUE)))
                                                                         ),

                                                                         # Information about number of tracks chosen:
                                                                         span(tags$b(textOutput("txtHubsize")), style="color:blue"),
                                                                         "Note: Some sample by assay combinations may either have both imputed and observed data or neither.",

                                                                         hr(),
                                                                         downloadButton('downSampleWUSTL', label="Download WUSTL TrackHub"),
                                                                         downloadButton('downSampleUCSC', label="Download UCSC TrackHub"),
                                                                         downloadButton('downSampleExport', label="Download File List"),
                                                                         h3("Usage Instructions:"),
                                                                         h4("WUSTL Epigenome Browser (Legacy):"),
                                                                         '1. Select tracks and download WUSTL TrackHub (json formatted)',
                                                                         br(),
                                                                         '2. Go to the ', tags$a(href="http://epigenomegateway.wustl.edu/legacy/?genome=hg19", "WUSTL Browser (legacy)"), 
                                                                         br(),
                                                                         '3. Select the Tracks > Custom Tracks > Add new tracks > Datahub by upload > Upload File',
                                                                         br(), br(),
                                                                         tags$b("NOTE: We do not recommend loading more than 100 tracks in the same WUSTL legacy track hub."),
                                                                         br(),
                                                                         "The WUSTL legacy client does not allow you to select specific tracks from a hub.",
                                                                         hr(),

                                                                         h4("WUSTL Epigenome Browser (new):"),
                                                                         '1. Select tracks and download WUSTL TrackHub (json formatted)',
                                                                         br(),
                                                                         "2. Go to the ", tags$a(href="http://epigenomegateway.wustl.edu/browser/?genome=hg19", "WUSTL Browser (new)"),
                                                                         br(),
                                                                         '3. Select the Tracks > Remote Tracks > Add Remote Data Hub > Upload File',
                                                                         br(), br(),
                                                                         "The new WUSTL browser does allow you to select specific tracks from a hub, so your hub can be much larger than the legacy version, as long as you do not load everything in.",
                                                                         br(),
                                                                         tags$b("Do not load in the coordinate override or native tracks - these are for the legacy version."),
                                                                         br(),

                                                                         hr(),
                                                                         h4("UCSC Epigenome Browser:"),

                                                                         '1. Select tracks and download UCSC TrackHub (plain text, all in one file)',
                                                                         br(),
                                                                         "2. Host hub on a http or ftp server (or use GBiB to host UCSC locally)",
                                                                         br(), 
                                                                         "3. Edit the following link with your hub location, then use it to connect the hub to the UCSC browser:",
                                                                         br(),
                                                                         tags$ul("http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&hubUrl=https://www.yourserver.org/hubfile.txt"),
                                                                         br(),
                                                                         tags$b("We do not recommend loading more than 1000 tracks in the same UCSC track hub."),
                                                                         br(), 
                                                                         "You do not have to display all tracks on the UCSC trackhub, but it may take a while to load in all of the data for many tracks."
                                                                )
                                                                ),
                                          )
                            )
                            ),
                   tabPanel("Modules/Motifs", value='enhNV',
                            sidebarPanel(selectizeInput('selMotifGroups', 'View specific sample group subpanel:',
                                                        choices = groups.full, multiple = FALSE, selected="Heart", options = list(create = TRUE)),
                                         hr(),
                                         radioButtons("motifDataset", "Motif dataset:",
                                                      c("Archetypes" = "arch", "All Motifs" = "full"), selected='arch'),
                                         checkboxInput("motifImages", "Show motifs as logos in network", value=FALSE),
                                         sliderInput("motifFCFilter", "Filter enrichments with absolute log2FC above:",
                                                     min = 0, max = 3, value = 1.5, step = .1)
                            ),


                            mainPanel(tabsetPanel(type="tabs", id="modulesTabset",
                                                  tabPanel("Fig. 2 Subpanels", value='tabMFigure', 
                                                           HTML("<div style ='overflow:auto; width:1000px ' >"),
                                                           htmlOutput("motifExample"),
                                                           HTML("</div>")
                                                  ),
                                                  tabPanel("Motif-Module Network", value='tabMNet', 
                                                           tags$b("Note:"), " Select nodes to see edge and node descriptions below the network.",
                                                           fluidRow(style='height:600px', column(12, withSpinner(visNetworkOutput("motifNetwork", height="600px"), color='#212121'))),
                                                           tags$b("Network edges (for selected node):"),
                                                           withSpinner(DT::dataTableOutput("networkSummaryTable"), color='#212121'),
                                                           hr(),
                                                           tags$b("Node summary (samples in cluster/motif logo):"),
                                                           withSpinner(DT::dataTableOutput("nodeSummaryTable"), color='#212121')
                                                           ),
                                                  tabPanel("Motif Enr. Table", value='tabMTable',
                                                           withSpinner(DT::dataTableOutput("MotifEnrTable"), color='#212121'),
                                                           downloadButton('downMotifEnrTable', label="Download Table")
                                                  ),
                                                  # NOTE: Too slow, consider taking out:
                                                  tabPanel("Module Centers", value='tabMMod', 
                                                           tags$b("Note:"), "Interactive heatmap (size = 833 x 300) may be slow to load.",
                                                           withSpinner(plotlyOutput("centersHeatmap", height = "1000px"), color="#212121")),
                                                  tabPanel("Motif Heatmap", value='tabMMotifHeat', 
                                                           tags$b("Note:"), "Interactive heatmap may be slow to load.",
                                                           withSpinner(plotlyOutput("motifHeatmap", height = "1000px"), color="#212121"))

                                         )
                            )),
                   tabPanel("GWAS Enrichments", value='treeNV',
                            useShinyjs(),
                            sidebarPanel(selectizeInput('selGWASTree', 'Select single GWAS of interest (PubMedID - Trait):',
                                                        choices = tree.gwas, multiple = FALSE, options = list(create = TRUE)),
                                         fluidRow(
                                                  column(4, actionButton("prevButton", "Previous"), actionButton("nextButton", "Next")),
                                                  column(1, tags$b("adj. p-value:")), 
                                                  column(4, selectizeInput('selFDR', label=NULL,
                                                                           # 'Select FDR (only for enr. vis):', 
                                                                           choices = c('1%', '0.1%'),
                                                                           multiple = FALSE, selected = '1%', options = list(create = TRUE)))
                                                  ),
                                         hr(),
                                         selectizeInput('filtNode', 'Filter tables by enrichment:',
                                                        choices = NULL, multiple = TRUE, options = list(create = TRUE)),
                                         selectizeInput('filtSNP', 'Filter tables by SNP:',
                                                        choices = NULL, multiple = TRUE, options = list(create = TRUE)),
                                         hr(),
                                         sliderInput("GWASsnpENHdist", "Maximum distance between SNP and center of enhancer:",
                                                     min = 0, max = 2500, value = 2500, step = 10),
                                         radioButtons("selGWASPrioType", "Show prioritization:",
                                                      c("All" = "all", "Has linked gene" = "any", "Linked gene disagrees with nearest gene" = "disagree"),
                                                      selected='all'),
                                         hr(),
                                         selectizeInput('selGWASTreeMult', 'Select multiple GWAS of interest to plot side-by-side:',
                                                        choices = tree.gwas, multiple = TRUE, options = list(create = TRUE)),
                                         sliderInput("GWASTreeScale", "Image Scale (for side-by-side):",
                                                     min = 0.25, max = 2, value = 1, step = 0.05)
                            ),

                            mainPanel(tabsetPanel(type="tabs", id="treeTabset",
                                                  tabPanel("Tree Overview", value='tabTree', 
                                                           tags$b("GWAS Enrichment on ", 
                                                                  tags$a(href=paste0(resdir, 'gwas_resources/enhancer_sample_tree.pdf'), "epigenomes hierarchy:")),
                                                           htmlOutput("GWASTreePng"), br(), uiOutput('renderGWASOverview')),
                                                  tabPanel("Enrichments", value='tabEnr', br(), uiOutput('renderGWASOverview2'),
                                                           hr(), DT::dataTableOutput("GWASOvlStats")),
                                                  tabPanel("Enhancers", value='tabEnh', br(), tags$b("Tissue-specific enhancers near GWAS lead SNPs:"),
                                                           "All enhancers within 2.5kb of a GWAS lead SNP that are also active ",
                                                           "in one of the top-enriched tree nodes in the GWAS.",
                                                           hr(), withSpinner(DT::dataTableOutput("GWASEnhStats"), color='#212121'),
                                                           downloadButton('downGWASEnh', label="Download Enhancers")),
                                                  tabPanel("Links", value='tabLinks', br(), tags$b("Table of gene-enhancer links:"), 
                                                           "Gene-enhancer links in the GWAS loci (SNPs +/- 1Mb), ",
                                                           "reported for the top-enriched sample groups in the GWAS.",
                                                           hr(), withSpinner(DT::dataTableOutput("GWASLinksStats"), color='#212121'),
                                                           downloadButton('downGWASLinks', label="Download Links")
                                                           ),
                                                  tabPanel("Enr. Heatmap", value='tabHeatmap',
                                                           
                                                           htmlOutput("GWASEnhancerSNP")),
                                                  tabPanel("Locus Vis.", value='tabLocusVis',
                                                           # Selection + navigation:
                                                           selectizeInput('selGWASLocus', 'Select locus:', choices = NULL, multiple=FALSE, options=list(create = TRUE)),
                                                           actionButton("prevLocusButton", "Previous"), 
                                                           actionButton("nextLocusButton", "Next"),
                                                           # For resizing GWAS locus:
                                                           tags$head(tags$style(type="text/css", "#GWASLocusPng img {max-width: 100%; width: 100%; height: auto}")),
                                                           tags$head(tags$script(src = "jquery.elevatezoom.min.js")),
                                                           tags$head(tags$script('Shiny.addCustomMessageHandler("zoomon",
                                                                                 function(message) {
                                                                                     var zoomConfig = {cursor: "crosshair", scrollZoom : true, 
                                                                                     zoomWindowPosition: 5, zoomWindowHeight:400, zoomWindowWidth:800}; 
                                                                                     $("#GWASLocusPng img").elevateZoom(zoomConfig); });')),
                                                           tags$head(tags$script('Shiny.addCustomMessageHandler("zoomoff",
                                                                                 function(message){
                                                                                     var image = $("#GWASLocusPng img");
                                                                                     // Remove old instance of EZ
                                                                                     $(".zoomContainer").remove();
                                                                                     image.removeData("elevateZoom");
                                                                                     image.removeData("zoomImage");
                                                                                 });')),

                                                           div(id='locusDiv', uiOutput("GWASLocusPng")),
                                                           htmlOutput("renderLocusOverview"),
                                                           tags$b("Click to enable/disable zoom on locus, scroll to change zoom size.")
                                                           ), 
                                                  tabPanel("Side-by-side", 
                                                           tags$b("Select multiple GWAS on sidebar to visualize their enrichments side by side:"),
                                                           uiOutput("GWASTreePngMultiple"))
                                                  ))
                            ),
                   tabPanel("About/Download", 
                            includeMarkdown("download.md"))
        )
)

