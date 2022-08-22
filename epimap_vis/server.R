#!/usr/bin/R 
# Server side for tree visualization:
library(tidyr)
library(shiny)
library(shinyTree)
library(stringr)
library(RColorBrewer)
library(visNetwork)
library(igraph)
library(heatmaply)
# library(plotly)
options(scipen=45)

source('auxiliary_plot_gwas.R')
source('auxiliary_trackhubs_functions.R')

# -------------
# Load metadata
# -------------
load('data/metadata_app.Rda') # Overall epi metadata
load('data/gwascatalog_may03_2019_noquotes_summary.Rda') # GWAS summaries
# For tables:
meta$ct = gsub("_"," ",meta$ct)
meta$infostr = paste0(meta$id, ': ', meta$info)
# For plots:
metamat = cbind('lifestage'=meta$lifestage, 'sex'=meta$sex,
                'type'=meta$type, 'project'=meta$Project,
                'group'=as.character(meta$GROUP))
meta$type[is.na(meta$type)] = 'ENCODE'
rownames(metamat) = meta$id

# ---------------------------
# Reduce to kept 859 samples:
# ---------------------------
cellorder = scan('data/833_samples.txt',quiet=TRUE, 'c', sep="\n")
co859 = scan('data/859_samples.txt',quiet=TRUE, 'c', sep="\n")
tree.gwas = scan('data/tree_gwas_list.txt', quiet=TRUE, 'c', sep="\n")
tree.gwas1p = scan('data/tree_gwas1p_list.txt', quiet=TRUE, 'c', sep="\n")
meta = meta[meta$id %in% co859,]

# Make tierdf:
all.assays = scan('data/all_assays.txt', quiet=TRUE, 'c', sep="\n")
mainmarks = c('DNase-seq', 'H3K27ac','H3K27me3','H3K36me3','H3K4me1','H3K4me3','H3K9me3')
tier2marks = c("ATAC-seq", "H2AFZ", "H3K4me2", "H3K79me2", "H3K9ac", "H4K20me1")
tier3marks = c("POLR2A", "CTCF", "EP300", "RAD21", "SMC3")
tierdf = data.frame(mark = all.assays, tier = 't4')
tierdf$tier[tierdf$mark %in% mainmarks] = 't1'
tierdf$tier[tierdf$mark %in% tier2marks] = 't2'
tierdf$tier[tierdf$mark %in% tier3marks] = 't3'

# For metadata tables:
mcols = c('GROUP','infoline','ct','lifestage','age','sex','type','Project','Donor')
rename.mcols = c('Group','Short name','Full name','Lifestage','Age','Sex','Type','Project','Donor')
# For colors/ordering:
odf = read.delim('data/group_ordering.tsv', header=T, stringsAsFactors=F)
clvs = c('Tissue', 'Primary Cell', 'Immune', 'Neural', 'Stem-like', 'Other')
odf$category = factor(odf$category, levels=clvs)
rdcol = unique(meta[,c('GROUP','COLOR')])
rdcol = rdcol[order(rdcol$GROUP),]
colset = as.character(rdcol$COLOR)
odf$alt = sub(" \\& ","/", odf$GROUP)

# Mapping from colors <-> groups
col2group = odf$GROUP
col2alt = odf$alt
group2col = odf$COLOR
names(col2group) = odf$COLOR
names(col2alt) = odf$COLOR
names(group2col) = odf$GROUP

# URLs:
epidir = 'https://personal.broadinstitute.org/cboix/epimap/'

# For loading RDa with no duplicates:
rda2list <- function(file) {
    e <- new.env()
    load(file, envir = e)
    as.list(e)
}


# Reactive values for keeping track of images:
imgTrack <- reactiveValues(newLocusImg=0)

shinyServer(function(input, output, session) {
  log <- c(paste0(Sys.time(), ": Interact with the tree to see the logs here..."))
  
  output$tree <- renderTree({
      if (input$selTree == 'groups'){
          # Tree from groups:
          mt2 = meta[meta$id %in% cellorder,]
          groups = sort(as.character(unique(mt2$GROUP)))
          grouptree = lapply(1:length(groups), function(i){
                                 x = groups[i]
                                 groupstr = tolower(gsub("[& -.]","",x))
                                 y = mt2$infostr[mt2$GROUP == x]
                                 ylist = lapply(1:length(y), function(j){
                                                    structure(y[j], stid=j, stclass='sample',
                                                              sttype=groupstr)})
                                 names(ylist) = y
                                 structure(ylist, stid=i, stclass='group', sttype=groupstr)})
          names(grouptree) = groups
          structure(grouptree, stopened=TRUE)
      } else {
          # Tree from enh:
          nlist = makeEnhTree()
          structure(nlist, stopened=TRUE)
      }
  })


  makeEnhTree <- reactive({
      load('data/enhancer_tree_declist.Rda')
      # Load in subsets:
      bdf = read.delim('data/prelim_e4_bkpts_loc.tsv', header=F)
      names(bdf) = c('node','name')
      bdf$uqname = make.unique(bdf$name, sep="_")
      # Collapse parents to create a collapsed tree:
      parents = declist$parent
      bind = !(parents %in% c(bdf$node, 1))
      while (sum(bind) > 0){
          parents[bind] = parents[parents[bind]]
          bind = !(parents %in% c(bdf$node, 1))
      }

      # Leaves:
      lind = which(declist$isleaf == 1)
      lmap = sapply(lind, function(x){ as.numeric(sub("_.*", "", declist$dec[[x]])) })
      # Show which leaves under what nodes:
      ldf = data.frame(ind=lmap,
                       lind=lind,
                       node=parents[lind],
                       id=declist$leaves[lmap],
                       nam=sapply(lind, function(x){declist$dec[[x]]}))
      ldf = ldf[order(ldf$ind),]

      # Create the collapsed tree:
      ldf$level = (parents[ldf$node] != 1) + (parents[parents[ldf$node]] != 1) + (ldf$node != 1)
      bdf$level = (parents[bdf$node] != 1) + (parents[parents[bdf$node]] != 1) + (bdf$node != 1)
      ldf = merge(ldf, bdf[,c('node','uqname')], all.x=TRUE)
      ldf$node[is.na(ldf$uqname)] = 1
      ldf$uqname[is.na(ldf$uqname)] = 'Misc.'
      ldf = ldf[order(ldf$node),]

      resdf = ldf
      # Build:
      nlist = list()
      bdf = bdf[order(bdf$level, decreasing=T),]
      for (i in 1:length(bdf$node)){
          node = bdf$node[i]
          nam = gsub("_", " ", bdf[bdf$node == node, 'uqname'])
          nids = ldf[ldf$node == node, 'id']
          y = meta$infostr[meta$id %in% nids]
          ygroup = meta$GROUP[meta$id %in% nids]
          groupstr = tolower(gsub("[& -.]","",ygroup))
          ids = sapply(nids, function(x){which(meta$id ==x)})
          names(ids) = NULL
          print(y)
          ylist = lapply(1:length(y), function(j){
                             structure(y[j], stid=ids[j] + nrow(bdf), stclass='sample',
                                       sttype=groupstr[j])})
          names(ylist) = y
          if (bdf$level[bdf$node == node] == 1){
              # Add the subset:
              subind = which(parents[bdf$node] == node)
              for (ind in rev(subind)){
                  snam = gsub("_", " ", bdf[ind, 'uqname'])
                  ylist2 = c(list(nlist[[snam]]), ylist)
                  names(ylist2)[1] = snam
                  print(ylist2)
                  ylist = ylist2
                  nlist[[snam]] = NULL
              }
          }
          structure(ylist, stid=i, stclass='group', sttype=groupstr[1])
          nlist[[nam]] = structure(ylist, stid=i, stclass='group', sttype=groupstr[1])
      }

      return(nlist)
  })


  getNodes <- reactive({
      tree <- input$tree
      text <- input$textNodes
      textCT <- input$textNodesCT
      nodes = NULL
      if (!is.null(tree)){ nodes = unlist(get_selected(input$tree)) }
      if (!is.null(text)){ nodes = sort(unique(c(nodes, text))) }
      if (!is.null(textCT)){ nodes = sort(unique(c(nodes, textCT))) }
      return(nodes)
  })

  output$selTxt <- renderPrint({
      nodes = getNodes()
      if (!is.null(nodes)){ 
          ids = unique(sub(": .*", "", nodes))
          paste0("https://epigenome.wustl.edu/epimap/data/imputed/impute_", ids, "_H3K27ac_DNase-seq.bigwig")
      }
  })

  getSampleTable <- reactive({
      nodes = getNodes()
      if (is.null(nodes)){ nodes = meta$infostr }
      ids = unique(sub(": .*", "", nodes))
      tab = meta[meta$id %in% ids,mcols]
      colnames(tab) = rename.mcols 
      tab
  })

  output$selTable = DT::renderDataTable({
      dt <- DT::datatable(data=getSampleTable(), options=list(pageLength=10))
      dt <- dt %>% formatStyle('Group', color='white')
      for (i in 1:nrow(odf)){
          dt <- dt %>% formatStyle('Group', backgroundColor=styleEqual(odf$GROUP[i], odf$COLOR[i]))
      }
      dt
  }, rownames=FALSE)

  output$downSampleTable <- downloadHandler(filename = 'metadata_subset.tsv',
                                            content = function(file){
                                                write.table(getSampleTable(), file, row.names=F, col.names=T, quote=F, sep="\t") 
                                            })

  # -----------------
  # Preset trackhubs:
  # -----------------
  getPresetURLs <- reactive({
      group = input$selTrackHubGroups
      groupstr = tolower(gsub("__", "", str_replace_all(group, "[^[:alnum:]]", "_")))
      # WUSTLurlprefix = "http://epigenomegateway.wustl.edu/legacy/?genome=hg19&datahub="
      WUSTLurlprefix="http://epigenomegateway.wustl.edu/browser/?genome=hg19&hub="
      UCSCurlprefix = "http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&hubUrl="
      hubprefix = paste0(epidir, 'trackhubs/pergroup/')
      WUSTLurlout = paste0(WUSTLurlprefix, hubprefix, 'wustlhub.fulldataset.pergroup.', groupstr, '.json')
      UCSCurlout = paste0(UCSCurlprefix, hubprefix, 'trackHub_bwtracks_pergroup_', groupstr, '.txt')
      return(list(WUSTL=WUSTLurlout, UCSC=UCSCurlout))
  })

  output$presetHubGroup <- renderUI({
      ll = getPresetURLs()
      urlWUSTL <- a("WUSTL", href=ll$WUSTL)
      urlUCSC <- a("UCSC", href=ll$UCSC)
      tagList(urlWUSTL, ' / ', urlUCSC)
  })

  # -----------------------------------------
  # Load the marks to export (for trackhubs):
  # -----------------------------------------
  getTrackList <- reactive({
      tdf = read.delim('data/all_released_tracks.tsv', header=F)
      names(tdf) = c('id','mark','prefix')
      rownames(meta) = meta$id
      tdf = merge(tdf, meta, all.x=TRUE)
      bwdir = 'https://epigenome.wustl.edu/epimap/data/'
      tdf$dir = paste0(bwdir, 'observed/')
      tdf$dir[grep('^impute', tdf$prefix)] = paste0(bwdir, 'imputed/')
      tdf$dataset = 'Observed'
      tdf$dataset[grep('^impute', tdf$prefix)] = 'Imputed'
      return(tdf)
  })

  getMarks <- reactive({
      marks = NULL
      if (length(input$selMarkTiers) > 0){
          # Tiers by data.frame:
          tiermarks = tierdf$mark[tierdf$tier %in% input$selMarkTiers]
          marks = c(marks, tiermarks) }
      if (length(input$selMarks) > 0){
          marks = c(marks, input$selMarks) }
      marks = sort(unique(marks))
      if (is.null(marks)){ marks = mainmarks }
      return(marks)
  })

  subsetTrackList <- reactive({
      # Get the subset of the tracklist + tier 1 marks only
      nodes = getNodes()
      ids = unique(sub(": .*", "", nodes))
      tdf = getTrackList()
      sdf = tdf[tdf$id %in% ids,]
      # Subset by marks:
      sdf = sdf[sdf$mark %in% getMarks(),]
      # Subset datatype:
      if (input$selIO == "imp"){
          sdf = sdf[sdf$dataset == "Imputed",]
      } else if (input$selIO == "obs"){
          sdf = sdf[sdf$dataset == "Observed",]
      }
      sdf
  })

  makeUCSCHub <- reactive({
      sdf = subsetTrackList()
      ucschub = make.ucsc(sdf)
      return(ucschub)
  })

  makeWUSTLHub <- reactive({
      sdf = subsetTrackList()
      jsonhub = make.legacy.wustl(sdf)
      return(jsonhub)
  })

  makeTracklist <- reactive({
      sdf = subsetTrackList()
      sdf$file = paste0(sdf$dir, sdf$prefix, 'wig.gz')
      sdf = sdf[,c('id','mark','file')]
      print(sdf)
      return(sdf)
  })


  output$txtHubsize = renderText({
      sdf = subsetTrackList()
      paste("Currently have", nrow(sdf), "tracks,",
            "covering", length(unique(sdf$id)), "samples",
            "across", length(unique(sdf$mark)), "assays")
  })

  output$downSampleWUSTL <- downloadHandler(filename = 'wustl_hub.config.json',
                                            content = function(file){
                                                write(makeWUSTLHub(), file) 
                                            })

  output$downSampleUCSC <- downloadHandler(filename = 'ucsc_trackHub.txt',
                                           content = function(file){
                                               write(makeUCSCHub(), file) 
                                           })

  output$downSampleExport <- downloadHandler(filename = 'track_filelist.txt',
                                             content = function(file){
                                                 write.table(makeTracklist(), file, row.names=F, col.names=T, quote=F, sep="\t") 
                                           })


  # ----------------------
  # Modules short example:
  # ----------------------
  mdir = "https://personal.broadinstitute.org/cboix/motif_examples/"

  output$motifExample <- renderText({
      mdir = "https://personal.broadinstitute.org/cboix/motif_examples/"
      if (length(input$selMotifGroups) == 1){
          groupstr = tolower(gsub("__", "", str_replace_all(input$selMotifGroups, "[^[:alnum:]]", "_")))
          imgsrc = paste0(mdir, 'clusters_motif_updated_groupdist_snapshot_group_', groupstr, '.png')
          c('<img src="',imgsrc,'" width=1000px >')
      }
  })

  getCentersData <- reactive({ 
      ll = rda2list('data/centers_enrichmat_data.Rda')
      return(ll)
  })

  # interactive heatmap prep
  setupCentersData <- reactive({
      ll = getCentersData()
      ll2 = getMotifData()
      colb <- colorRampPalette(brewer.pal(n=9,name="Blues"))(100)
      mat = ll$centers[rev(ll$rngroup), ll$clsord]
      rownames(mat) = paste0(rownames(mat),": ", meta[rownames(mat),'infoline'])
      
      rfac = meta[rev(ll$rngroup),'COLOR']
      cfac = ll2$clsmaxcols[ll$clsord]
      rfac = factor(col2group[rfac])
      cfac = factor(col2group[cfac])
      heatmaply(mat, colors=colb, limits=c(0,1),
                row_side_colors=rfac, col_side_colors=cfac,
                row_side_palette=group2col, col_side_palette=group2col,
                hide_colorbar=TRUE,
                side_color_colorbar_len=0,
                Rowv=FALSE, Colv=FALSE, srtCol = 45) %>% 
          layout(margin = list(l = 100, b = 200))
  })

  # interactive heatmap prep
  setupMotifHeatmapData <- reactive({
      ll = getCentersData()
      ll2 = getMotifData()
      colrb <- colorRampPalette(brewer.pal(n=7,name="RdBu"))(100)
      if (input$motifDataset == 'arch'){
          mat = t(ll2$mmat[ll$clsord, rev(ll2$motnam)])
      } else {
          mat = t(ll2$full.mmat[ll$clsord, rev(ll2$full.motnam)])
      }
      zlim = 1.5
      mat[mat > zlim] = zlim
      mat[mat < -zlim] = -zlim
      cfac = ll2$clsmaxcols[ll$clsord]
      cfac = factor(col2group[cfac])
      heatmaply(mat, colors=rev(colrb), limits=c(-zlim, zlim),
                col_side_colors=cfac, col_side_palette=group2col,
                side_color_colorbar_len=0,
                Rowv=FALSE, Colv=FALSE, srtCol = 45) %>% 
          layout(margin = list(l = 100, b = 200))
  })



  output$centersHeatmap <- renderPlotly({
      setupCentersData()
  })

  output$motifHeatmap <- renderPlotly({
      setupMotifHeatmapData()
  })

  getMotifData <- reactive({
      ll = rda2list('data/motif_modules_network_data.Rda')
      ll$rev.mmap = names(ll$mmap)
      ll$rev.full.mmap = names(ll$full.mmap)
      names(ll$rev.mmap) = ll$mmap
      names(ll$rev.full.mmap) = ll$full.mmap
      names(ll$rev.mmap) = sub('\\+', '.', names(ll$rev.mmap))
      names(ll$rev.full.mmap) = sub('\\+', '.', names(ll$rev.full.mmap))
      return(ll)
  })

  makeMotifReducedDF <- reactive({
      ll = getMotifData()
      mdf = data.frame(ll$mmat)
      mdf$cls = rownames(ll$mmat)
      mdf = gather(mdf, motif, log2FC, -cls)
      mdf = mdf[mdf$log2FC != 0,]
      mdf$log2FC = round(mdf$log2FC, 3)
      mdf$full.name = ll$rev.mmap[mdf$motif]
      mdf$cls.group = col2alt[ll$clsmaxcols[mdf$cls]]
      return(mdf)
  })

  makeMotifFullDF <- reactive({
      ll = getMotifData()
      mdf = data.frame(ll$full.mmat)
      mdf$cls = rownames(ll$full.mmat)
      mdf = gather(mdf, motif, log2FC, -cls)
      mdf = mdf[mdf$log2FC != 0,]
      mdf$log2FC = round(mdf$log2FC, 3)
      mdf$full.name = ll$rev.full.mmap[mdf$motif]
      mdf$cls.group = col2alt[ll$clsmaxcols[mdf$cls]]
      return(mdf)
  })

  # Potentially: 
  # Enable select specific groups of samples
  # or any cluster with X sample in it?
  subsetMotifDF <- reactive({
      if (input$motifDataset == "arch"){
          mdf = makeMotifReducedDF()
      } else {
          mdf = makeMotifFullDF()
      }
      if (input$motifFCFilter != 0){
          mdf = mdf[abs(mdf$log2FC) >= input$motifFCFilter,]
      }
      logodir = "https://personal.broadinstitute.org/cboix/motif_data_epimap/img/"
      mdf$motif.logo = paste0('<img src="', logodir, 'logo_', mdf$full.name, '.png" height="20"></img>')
      return(mdf)
  })

  output$MotifEnrTable = DT::renderDataTable({
      df = subsetMotifDF()
      dt = DT::datatable(df, escape = FALSE, rownames=FALSE) # To render images:
      print(odf$GROUP)
      # Colors:
      brks <- c(-5,seq(-2.5, 2.5, by=.25),5)
      clrs <- colorRampPalette(brewer.pal(n=9,name="RdBu"))(length(brks) + 1)
      dt <- dt %>% formatStyle('log2FC', color=styleInterval(c(-1.5, 1.5), c('white','black','white')))
      dt = dt %>% formatStyle('log2FC', backgroundColor=styleInterval(brks, rev(clrs)))
      dt <- dt %>% formatStyle('cls', 'cls.group', color=styleEqual(odf$alt, ifelse(odf$alt %in% c('Multiple', 'Eye','Lung'),'grey20', 'white')))
      dt <- dt %>% formatStyle('cls', 'cls.group', backgroundColor=styleEqual(odf$alt, odf$COLOR))
      dt
  })

  output$downMotifEnrTable <- downloadHandler(filename = 'motif_enrichments.tsv',
                                              content = function(file){
                                                  write.table(subsetMotifDF(), file, row.names=F, col.names=T, quote=F, sep="\t") 
                                              })

  makeMotifNetwork <- reactive({
      ll = getMotifData()
      # Process relative to options:
      if (input$motifDataset == 'arch'){
          imat = ll$mmat
      } else {
          imat = ll$full.mmat
      }
      imat[abs(imat) < input$motifFCFilter] = 0
      imat = imat
      color.edges = TRUE
      color.max = FALSE

      # Reduce to containing motif:
      imat = imat[apply(abs(imat) > 0, 1, sum) > 0,]
      imat = imat[,apply(abs(imat) > 0, 2, sum) > 0]
      dim(imat)

      # Should move all of this into a function in aux. funcs:
      # Make and modify graph:
      net <- graph_from_incidence_matrix(imat, weighted=TRUE)
      ewt = E(net)$weight
      E(net)$weight = 1
      V(net)$color <- c("steel blue", "orange")[V(net)$type+1]
      V(net)$shape <- c("square", "circle")[V(net)$type+1]
      V(net)$frame.color = NA
      # Labels + cols:
      lbnames = attributes(V(net))$names
      lind = lbnames %in% rownames(imat)
      if (color.max){
          V(net)$color[lind] = ll$clsmaxcols[lbnames[lind]]
      } else {
          V(net)$color[lind] = ll$clscols[lbnames[lind]]
      }
      V(net)$label[lind] = ""
      V(net)$label[!lind] = lbnames[!lind]
      V(net)$label.cex=.4
      V(net)$label.font=2
      V(net)$label.family=2

      E(net)$width <- abs(ewt) / 4
      ewt[ewt > 3] = 3
      ewt[ewt < -3] = -3

      colrb <- colorRampPalette(brewer.pal(n=7,name="RdBu"))(100)
      col_fun = function(x, pal=rev(colrb[10:90])){
          bin <- cut(x, seq(-3, 3, length.out=length(pal)), include.lowest=T) 
          pal[bin] }
      if (color.edges){ E(net)$color = col_fun(ewt) }
      set.seed(1)

      return(net=net)
  })

  output$motifNetwork <- renderVisNetwork({
      net = makeMotifNetwork()
      ll = getMotifData()
      data <- toVisNetworkData(net)
      # If use motif logo images:
      if (input$motifImages){
          logodir = "https://personal.broadinstitute.org/cboix/motif_data_epimap/img/"
          nodnam = data$nodes$id
          nind = grep("^c[0-9]*$",nodnam, invert=TRUE)
          data$nodes$image = ""
          if (input$motifDataset == 'arch'){
              fnam = ll$rev.mmap[nodnam[nind]]
          } else {
              fnam = ll$rev.full.mmap[nodnam[nind]]
          }
          data$nodes$shape[nind] = "image"
          data$nodes$image[nind] = paste0(logodir, 'logo_', fnam, '.png')
      }
      # logodir, 'logo_', df$full.name, '.png'
      visNetwork(nodes = data$nodes, edges = data$edges, width='1000px', height='1000px') %>%
          visIgraphLayout(layout='layout_with_fr') %>%
          visInteraction(hover = TRUE) %>%
          visOptions(highlightNearest = TRUE,
                     nodesIdSelection = list(enabled = TRUE)) %>%
          visPhysics(stabilization = FALSE) %>%
          visEdges(smooth = FALSE) %>%
          visEvents(hoverNode = "function(nodes) { Shiny.onInputChange('current_node_id', nodes); ;}") %>%
          visInteraction(navigationButtons = TRUE)
  })

  output$networkSummaryTable = DT::renderDataTable({
      sel.node = input$motifNetwork_selected
      sel.node = sub("\\+", ".", sel.node) # For GATA+TAL1 and similar.
      df = subsetMotifDF()
      node.type = "motif"
      if (sel.node != ""){
          # Select node or cluster:
          if (length(grep("^c[0-9]*$",sel.node)) > 0){
              node.type = "cluster" 
              df = df[df$cls %in% sel.node,]
              df = df[order(df$log2FC, decreasing=T),]
              dt = DT::datatable(df, escape = FALSE, rownames=FALSE)
          } else {
              ll = getMotifData()
              df = df[df$motif %in% sel.node,]
              df = df[order(df$log2FC, decreasing=T),]
              dt = DT::datatable(df, escape = FALSE, rownames=FALSE)
          }
          # Colors:
          brks <- c(-5,seq(-2.5, 2.5, by=.25),5)
          clrs <- colorRampPalette(brewer.pal(n=9,name="RdBu"))(length(brks) + 1)
          dt <- dt %>% formatStyle('log2FC', color=styleInterval(c(-1.5, 1.5), c('white','black','white')))
          dt <- dt %>% formatStyle('log2FC', backgroundColor=styleInterval(brks, rev(clrs)))
          dt <- dt %>% formatStyle('cls', 'cls.group', color=styleEqual(odf$alt, ifelse(odf$alt %in% c('Multiple', 'Eye','Lung'),'grey20', 'white')))
          dt <- dt %>% formatStyle('cls', 'cls.group', backgroundColor=styleEqual(odf$alt, odf$COLOR))
      } else {
          df = data.frame()
          dt = DT::datatable(df, escape = FALSE, rownames=FALSE)
      }
      dt
  })

  output$nodeSummaryTable = DT::renderDataTable({
      sel.node = input$motifNetwork_selected
      sel.node = sub("\\+", ".", sel.node) # For GATA+TAL1 and similar.
      node.type = "motif"
      if (length(grep("^c[0-9]*$",sel.node)) > 0){
          node.type = "cluster" 
          ll = getCentersData()
          xc = ll$centers[, sel.node]
          xind = (xc >=.25)
          df = data.frame(cls=sel.node, id=rownames(ll$centers)[xind], inclusion=round(xc[xind],2))
          df$name = meta[df$id, 'infoline']
          df$group = sub(" \\& ", "/", meta[df$id, 'GROUP'])
          df$type = meta[df$id,'type']
          df$lifestage = meta[df$id,'lifestage']
          df = df[order(df$inclusion, decreasing=T),]
          # Color/format:
          dt = DT::datatable(df, escape = FALSE, rownames=FALSE) # To render images:
          dt <- dt %>% formatStyle('name', 'group', color=styleEqual(odf$alt, ifelse(odf$alt == 'Multiple','grey20', 'white')))
          dt <- dt %>% formatStyle('name', 'group', backgroundColor=styleEqual(odf$alt, odf$COLOR))
          dt <- dt %>% formatStyle('inclusion', color='white')
          dt = dt %>% formatStyle('inclusion', background = styleColorBar(c(0,1), 'darkblue'),
                                  # backgroundSize = '98% 88%',
                                  backgroundSize = '98% 75%',
                                  backgroundRepeat = 'no-repeat',
                                  backgroundPosition = 'center')
      } else if (sel.node != ""){
          df = subsetMotifDF()
          df = df[df$motif %in% sel.node,]
          df = df[order(df$log2FC, decreasing=T),]
          df = df[1,c('motif','full.name','motif.logo')]
          dt = DT::datatable(df, escape = FALSE, rownames=FALSE)
      } else {
          df = data.frame()
          dt = DT::datatable(df, escape = FALSE, rownames=FALSE)
      }
      dt
  })



  # ---------------------------------
  # Sample tree and gwas enrichments:
  # ---------------------------------
  getTreePng <- reactive({
      trait = sub(".* - ", "", input$selGWASTree)
      pmid = sub(" - .*", "", input$selGWASTree)
      traitstr = gsub("'","_", gsub(" ", "_", tolower(trait)))
      traitstrnoparen = gsub("/", "_", gsub("\\)" ,"", gsub("\\(","",traitstr)))
      pnghead = paste0(epidir, 'gwas_smallfigures/')
      pngpref = 'enhancers_e2500_adj_onecut_cons_parent_gwas_small_'
      if (input$selFDR == '1%'){ pcutstr = '_1pct' } else { pcutstr = '_pt1pct' }
      pngurl = paste0(pnghead, pngpref, traitstrnoparen, '_', pmid, pcutstr, '.png')
      return(pngurl)
  })

  output$GWASTreePng <- renderText({
      imgsrc <- getTreePng()
      c('<img src="',imgsrc,'" width=600px >')
  })

  getEnhancerSNPPng <- reactive({
      ind = which(tree.gwas == input$selGWASTree)
      ipref = sprintf("%05d", ind)
      pnghead = paste0(epidir, 'gwas_smallfigures/')
      pngpref = 'enhancers_e2500_adj_onecut_snp_intersections_'
      pngurl = paste0(pnghead, pngpref, sprintf('page_%d.png', ind))
  })

  output$GWASEnhancerSNP <- renderText({
      imgsrc <- getEnhancerSNPPng()
      c('<img src="',imgsrc,'" width=500px >')
  })

  getGWASList <- reactive({
      if (input$selFDR == '1%'){ tree.gwas } else { tree.gwas1p }
  })

  observeEvent(input$selFDR, {
                   tree.list = getGWASList()
                   updateSelectizeInput(session, "selGWASTree", choices=tree.list)
                   updateSelectizeInput(session, "selGWASTreeMult", choices=tree.list)
  })

  # Navigate GWAS selection with buttons
  bnext <- observeEvent(input$nextButton, {
                            tree.list = getGWASList()
                            ind = which(tree.list == input$selGWASTree)
                            if (ind != length(tree.list)){ ind <- ind + 1 } 
                            updateSelectizeInput(session, "selGWASTree", selected=tree.list[ind]) })

  bprev <- observeEvent(input$prevButton, {
                            tree.list = getGWASList()
                            ind = which(tree.list == input$selGWASTree)
                            if (ind != 1){ ind <- ind - 1 } 
                            updateSelectizeInput(session, "selGWASTree", selected=tree.list[ind]) })

  getTreeMultiPNG <- reactive({
      tree.list = getGWASList()
      trait = sub(".* - ", "", input$selGWASTreeMult)
      pmid = sub(" - .*", "", input$selGWASTreeMult)
      traitstr = gsub("'","_", gsub(" ", "_", tolower(trait)))
      traitstrnoparen = gsub("/", "_", gsub("\\)" ,"", gsub("\\(","",traitstr)))
      pnghead = paste0(epidir, 'gwas_smallfigures/')
      pngpref = 'enhancers_e2500_adj_onecut_cons_parent_gwas_small_'
      if (input$selFDR == '1%'){ pcutstr = '_1pct' } else { pcutstr = '_pt1pct' }
      pngurl = paste0(pnghead, pngpref, traitstrnoparen, '_', pmid, pcutstr, '.png')
      return(pngurl)
  })

  getTreeMultiScale <- reactive({
      return(input$GWASTreeScale * 600)
  })

  output$GWASTreePngMultiple <- renderText({
      imgsrc <- getTreeMultiPNG()
      imgscale  <- getTreeMultiScale()
      paste('<img src="',imgsrc,'" width=', imgscale, 'px >')
  })

  # ----------------------------------------
  # Tables for the GWAS Tree prioritization:
  # ----------------------------------------
  # Overview of GWAS:
  # ind = which(gwas.summary.df$uid %in% c(input$selGWASTree, input$selGWASTreeMult))
  # gwas.summary.df[ind, c('pubMedID','trait','sampleSize','numSNPs','initSample','replSample','pubDate')]
  # Reactive changing data:
  # 'linking_data/all_snp_enrichment_labs_forfilt.Rda'

  # Dynamically update the slicing:
  getGWASFiltLabs <- reactive({
      load('data/all_snp_enrichment_labs_forfilt.Rda')
      return(list(snp=snplist, enr=enrlist, groups=linklist))
  })

  selectGWASFiltLabs <- reactive({
      ll = getGWASFiltLabs()
      return(list(snp = ll$snp[[input$selGWASTree]], enr = ll$enr[[input$selGWASTree]], groups=ll$groups[[input$selGWASTree]]))
  })

  observeEvent(input$selGWASTree, {
                   ll = selectGWASFiltLabs()
                   updateSelectizeInput(session, "filtNode", choices=ll$enr)
                   updateSelectizeInput(session, "filtSNP", choices=ll$snp)
                   updateSelectizeInput(session, "selGWASLocus", choices=ll$snp)
  })

  output$renderGWASOverview <- renderUI({
      df = gwas.summary.df[gwas.summary.df$uid == input$selGWASTree,]
      urlPM = strong(a(paste0("(PubMed: ", df$pubMedID[1],"):"), href=paste0("https://pubmed.ncbi.nlm.nih.gov/", df$pubMedID[1])))
      title = strong(paste0(df$trait[1], " GWAS"))
      desc = "Top-enriched epigenome tree nodes."
      ssamp = strong('Sample size:') 
      tagList(title, urlPM, desc, br(),
              ssamp, df$initSample[1], br(), 
              strong('Date:'), df$pubDate[1])
  })

  # Same as above, for another tab:
  output$renderGWASOverview2 <- renderUI({
      df = gwas.summary.df[gwas.summary.df$uid == input$selGWASTree,]
      urlPM = strong(a(paste0("(PubMed: ", df$pubMedID[1],"):"), href=paste0("https://pubmed.ncbi.nlm.nih.gov/", df$pubMedID[1])))
      title = strong(paste0(df$trait[1], " GWAS"))
      desc = "Top-enriched epigenome tree nodes and nominal enrichment p-values."
      ssamp = strong('Sample size:') 
      tagList(title, urlPM, desc, br(),
              ssamp, df$initSample[1], br(), 
              strong('Date:'), df$pubDate[1])
  })

  # Add also top enriched tissues, number of SNPs intersecting them, etc.
  getGWASOvlStats <- reactive({
      load('data/all_enrichment_summaries_alone.Rda')
      gw.enrdf
  })

  subsetGWASOvlStats <- reactive({
      df <- getGWASOvlStats()
      df <- df[df$uid %in% input$selGWASTree,]
      if (!is.null(input$filtNode)){
          df = df[df$enrlab %in% input$filtNode,]
      }
      df = df[order(df$node.rank),]
      if (input$selFDR == '0.1%'){ df = df[df$enr.p <= 10^(-22.09788),] } 
      df$enr.p = signif(df$enr.p, 2)
      df = df[, c('node.rank', 'node.name','enr.p', 'pubMedID','trait', 'node.group')]
      colnames(df) = c('enrRank','enrName','enr.pValue', 'pubMedID','trait', 'enrGroup')
      rownames(df) = NULL
      df
  })

  colorDT = function(df){
      dt <- DT::datatable(data=df)
      dt <- dt %>% formatStyle('enrName', 'enrGroup', color=styleEqual(odf$GROUP, ifelse(odf$GROUP == 'Multiple','grey20', 'white')))
      dt <- dt %>% formatStyle('enrName', 'enrGroup', backgroundColor=styleEqual(odf$GROUP, odf$COLOR))
      dt
  }

  output$GWASOvlStats = DT::renderDataTable({
      df = subsetGWASOvlStats()
      colorDT(df)
  }, rownames= FALSE)

  # All nearby enhancers:
  getGWASEnhStats <- reactive({
      load('data/all_snpxenhancer_intersections.Rda')
      snpenhdf
  })

  # Subset to the proper uid:
  subsetGWASEnhStats <- reactive({
      df = getGWASEnhStats()
      sdf = subsetGWASOvlStats()
      sdf = sdf[,c('enrName','enrGroup')]
      df = df[df$uid %in% input$selGWASTree,]
      if (input$selFDR == '0.1%'){ df = df[df$enr.p <= 10^(-22.09788),] } 
      df$enr.p = signif(df$enr.p,2)
      df = df[, c('chr','snpPos','p', 'start','end','dist','nearest','node.rank','node.name','enr.p', 'pubMedID','trait', 'enrlab', 'snplab')]
      colnames(df) = c('chr','snpPos','snp.pValue', 'enhStart','enhEnd','distToCenter','nearestGene','enrRank','enrName','enr.pValue', 'pubMedID','trait', 'enrlab','snplab')
      df = merge(df, sdf)
      df
  })

  # Filter down the data subsetted to the correct uid:
  filtGWASEnhStats <- reactive({
      df = subsetGWASEnhStats()
      if (!is.null(input$filtNode)){
          df = df[df$enrlab %in% input$filtNode,]
      }
      if (!is.null(input$filtSNP)){
          df = df[df$snplab %in% input$filtSNP,]
      }
      if (input$GWASsnpENHdist != 2500){
          df = df[df$distToCenter <= input$GWASsnpENHdist,]
      }
      df[, c('chr','snpPos','snp.pValue', 'enhStart','enhEnd','distToCenter','nearestGene','enrRank','enrName','enr.pValue', 'enrGroup')] #, 'pubMedID','trait')]
  })

  output$GWASEnhStats = DT::renderDataTable({
      df = filtGWASEnhStats()
      colorDT(df)
  }, rownames= FALSE)

  output$downGWASEnh <- downloadHandler(filename = 'active_enh_inloci.tsv',
                                            content = function(file){
                                                write.table(filtGWASEnhStats(), file, row.names=F, col.names=T, quote=F, sep="\t") 
                                            })

  getGWASLinksStats  <- reactive({
      load('data/all_gwas_SNP_links.Rda')
      gw.linksdf$score = round(gw.linksdf$score, 2)
      scols = c('chr','snpPos','p','dist','nearest','symbol','score','linkdist','node.rank','node.name', 'enr.p','uid', 'enrlab','snplab')
      final.scols = c('chr','snpPos','snp.pValue','distToCenter','nearestGene','linkedGene','linkScore','linkDist','enrRank','enrName','enr.pValue', 'uid', 'enrlab','snplab')
      gw.linksdf = gw.linksdf[,scols]
      colnames(gw.linksdf) = final.scols
      gw.linksdf
  })

  subsetGWASLinksStats <- reactive({
      df = getGWASLinksStats()
      df = df[df$uid == input$selGWASTree,]
      if (input$selGWASPrioType == 'disagree'){ 
          df = df[!is.na(df$linkedGene),]
          df = df[df$nearestGene != df$linkedGene,] 
      } else if (input$selGWASPrioType == 'any'){
          df = df[!is.na(df$linkedGene),]
      }
      df = df[df$distToCenter <= input$GWASsnpENHdist,]
      sdf = subsetGWASOvlStats()
      sdf = sdf[,c('enrName','enrGroup')]
      df = merge(df, sdf, all.x=TRUE)
      if (input$selFDR == '0.1%'){ df = df[df$enr.pValue <= 10^(-22.09788),] } 
      final.scols = c('chr','snpPos','snp.pValue','distToCenter','nearestGene','linkedGene','linkScore','linkDist',
                      'enrRank','enrName','enr.pValue', 'enrGroup', 'enrlab','snplab')
      unique(df[,final.scols])
  })

  filtGWASLinksStats <- reactive({
      df = subsetGWASLinksStats()
      if (!is.null(input$filtNode)){
          df = df[df$enrlab %in% input$filtNode,]
      }
      if (input$treeTabset == "tabLocusVis"){
          df = df[df$snplab %in% input$selGWASLocus,]
      } else {
          if (!is.null(input$filtSNP)){
              df = df[df$snplab %in% input$filtSNP,]
          }
      }
      if (input$GWASsnpENHdist != 2500){
          df = df[df$distToCenter <= input$GWASsnpENHdist,]
      }
      final.scols = c('chr','snpPos','snp.pValue','distToCenter','nearestGene','linkedGene','linkScore','linkDist','enrRank','enrName','enr.pValue', 'enrGroup')
      df = df[order(df$linkScore, decreasing=TRUE),] 
      df = df[order(df$enr.pValue),] 
      df = df[order(df$snp.pValue),] 
      df[, final.scols]
  })

  colorGWASLinksStats <- reactive({
      df = filtGWASLinksStats()
      dt = colorDT(df)
      # Also change distances:
      brks <- seq(0, 2500, by=500) 
      clrs <- colorRampPalette(brewer.pal(n=9,name="Reds"))(length(brks) + 1)
      brks2 <- seq(0, 1e6, by=50000) 
      clrs2 <- colorRampPalette(brewer.pal(n=9,name="Blues"))(length(brks2))
      brks2 = c(-rev(brks2), brks2[-1])
      clrs2 = c(clrs2, rev(clrs2))
      dt = dt %>% formatStyle('distToCenter', backgroundColor=styleInterval(brks, rev(clrs)))
      dt <- dt %>% formatStyle('linkDist', color=styleInterval(c(-1e6,-5e5,5e5,1e6), c('black','black','white','black', 'black')))
      dt = dt %>% formatStyle('linkDist', backgroundColor=styleInterval(brks2, clrs2))
      # Also change the scores to bars:
      dt = dt %>% formatStyle('linkScore', background = styleColorBar(c(0,1), 'lightblue'),
                              # backgroundSize = '98% 88%',
                              backgroundSize = '98% 75%',
                              backgroundRepeat = 'no-repeat',
                              backgroundPosition = 'center')
      dt
  })

  output$GWASLinksStats = DT::renderDataTable({
      colorGWASLinksStats()
  }, rownames= FALSE)


  output$downGWASLinks <- downloadHandler(filename = 'corrlinks_inloci.tsv',
                                            content = function(file){
                                                write.table(filtGWASLinksStats(), file, row.names=F, col.names=T, quote=F, sep="\t") 
                                            })

  # NOTE: Change to show the different table here, if we keep this section:
  output$GWASLinksStats2 = DT::renderDataTable({
      colorGWASLinksStats()
  }, rownames= FALSE)

  output$downGWASLinks2 <- downloadHandler(filename = 'corrlinks_inloci.tsv',
                                            content = function(file){
                                                write.table(filtGWASLinksStats(), file, row.names=F, col.names=T, quote=F, sep="\t") 
                                            })


  # -------------------------
  # GWAS Locus Visualization:
  # -------------------------
  # Navigate locus with buttons
  bnextlocus <- observeEvent(input$nextLocusButton, {
                                 # session$sendCustomMessage(type = 'zoomoff', message = list())
                                 ll = selectGWASFiltLabs()
                                 snpset = ll$snp
                                 ind = which(snpset == input$selGWASLocus)
                                 if (ind != length(snpset)){ ind <- ind + 1 } 
                                 updateSelectizeInput(session, "selGWASLocus", selected=snpset[ind]) 
                                            })

  bprevlocus <- observeEvent(input$prevLocusButton, {
                                 # Reset the zoom:
                                 # session$sendCustomMessage(type = 'zoomoff', message = list())
                                 ll = selectGWASFiltLabs()
                                 snpset = ll$snp
                                 ind = which(snpset == input$selGWASLocus)
                                 if (ind != 1){ ind <- ind - 1 } 
                                 updateSelectizeInput(session, "selGWASLocus", selected=snpset[ind])
                                            })


  getGWASLocusPNG <- reactive({
      # Directory for the loci:
      suid = input$selGWASTree
      trait = sub(".* - ", "", suid)
      pmid = sub(" - .*", "", suid)
      traitstr = gsub("'","_", gsub(" ", "_", tolower(trait)))
      traitstrnoparen = gsub("/", "_", gsub("\\)" ,"", gsub("\\(","",traitstr)))
      pnghead = paste0(epidir, 'gwas_smallfigures/gwas_locus_figures/')
      pngpref = paste0(traitstrnoparen, '_', pmid, '/gene_link_vis_')
      pngtail = '_w1000000_with_tx_grouponly.png'
      # Get locus:
      pnglocus = sub(":","_", gsub(",","", sub(" .*", "", input$selGWASLocus)))
      pngurl = paste0(pnghead, pngpref, pnglocus, pngtail)
      # chr6_43341689
      return(pngurl)
  })

  # output$GWASLocusPng <- renderText({
  #     imgsrc <- getGWASLocusPNG()
  #     c('<img src="',imgsrc,'" width=900px >')
  # })

  output$GWASLocusPng <- renderUI({
      imgsrc <- getGWASLocusPNG()
      session$sendCustomMessage(type = 'zoomoff', message = list()) 
      # Update reactive values, indicating a new img was loaded:
      imgTrack$newLocusImg = 1 
      # Image tags:
      img(src=imgsrc, "data-zoom-image"=imgsrc, width=1000)
  })

  switchZoom = function(){
      if (imgTrack$newLocusImg == 1){
          session$sendCustomMessage(type = 'zoomon', message = list())
          imgTrack$newLocusImg = 0
      } else {
          session$sendCustomMessage(type = 'zoomoff', message = list())
          imgTrack$newLocusImg = 1
      } 
  }

  shinyjs::onclick("GWASLocusPng", switchZoom())

  # Disable zoom if change tabs:
  observeEvent(input$treeTabset, {
                   if (input$treeTabset != "tabLocusVis"){
                       session$sendCustomMessage(type = 'zoomoff', message = list())
                       imgTrack$newLocusImg = 1
  }})

  # Disable zoom if on other pages:
  observeEvent(input$navBarSet, {
                   if (input$navBarSet != "treeNV"){
                       session$sendCustomMessage(type = 'zoomoff', message = list())
                       imgTrack$newLocusImg = 1
  }})

  output$renderLocusOverview <- renderText({
      ll = selectGWASFiltLabs()
      link.groups = ll$groups 
      link.cols = sapply(link.groups, function(x){ odf[odf$GROUP == x, 'COLOR']})
      txt.cols = ifelse(link.groups %in% c('Multiple','Eye','Lung','Other'), 'black','white')
      # 
      suid = input$selGWASTree
      trait = sub(".* - ", "", suid)
      pmid = sub(" - .*", "", suid)
      traitstr = gsub("'","_", gsub(" ", "_", tolower(trait)))
      traitstrnoparen = gsub("/", "_", gsub("\\)" ,"", gsub("\\(","",traitstr)))
      gwasdir = paste0(epidir, 'gwas_smallfigures/gwas_locus_figures/')
      traitdir = paste0(gwasdir, traitstrnoparen, '_', pmid, '/')
      # Datasets:
      traittsv = paste0(traitdir, traitstrnoparen, '_plotted_links_allsnps_w1000000.tsv.gz')
      traitrda = paste0(traitdir, traitstrnoparen, '_plotted_links_allsnps_w1000000.Rda')
      # df = gwas.summary.df[gwas.summary.df$uid == input$selGWASTree,]
      urlPM = strong(a(paste0("data repository"), href=paste0(traitdir)))

      # Create legend:
      fmt.groups = c()
      for (i in 1:length(link.groups)){
          fmt.groups = c(fmt.groups, paste0('<strong><span style="background-color:', link.cols[i],
                                            '; color:', txt.cols[i], '">', link.groups[i], "</span></strong>")) }
      fmt.groups = paste(fmt.groups, collapse=', ')
      title = strong("Locus overview for 1Mb around selected lead SNP:")
      linkdesc = tagList("Two types of correlation-based links are plotted:", strong("(1)"), "Links from one of the enhancers near a lead SNP in the enriched epigenomes.", strong("(2)"), "Any links in the locus present in at least half of the samples in one of the top sample groups (")
      taildesc = tagList( "Genes linked to an enhancer within 2.5kb of a GWAS lead SNP are highlighted and colored according to the sample group with the highest link score.")
      pltdesc = "Tracks show average H3K27ac signal of enhancers in locus, and red dashed lines indicate the TSSes of nearby genes."
      urldesc = "Link data and images for this GWAS are also available from our "
      paste0(tagList(title, pltdesc, linkdesc), fmt.groups, '). ',
             tagList(taildesc, urldesc, urlPM, "."))
  })

})

