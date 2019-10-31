#library(CGtools)

chrom_info <- list(human = paste("chr", c(1:22, "X"), sep=""), # hg18
                   mouse = paste("chr", c(1:19, "X"), sep=""), # mm9
                   fly   = c("2L", "2R", "3L", "3R", "4", "X"), # dm3
                   fugu  = "chrUn");

#chromsize_info <- list(human = chromSE_HSA_NCBI36, # hg18
#                       mouse = chromSE_MM_NCBIm37, # mm9
#                       fly   = chromSE_DMEL_R5,    #dm3
#                       fugu  = data.frame(seqname="chrUn", start=1, end=400509343));

id <- function() paste("WM", gsub("-", "", Sys.Date()), sep="")

#runmean <- function(x,k) {
#  if ((k %% 2) == 0) { k <- k+1; }
#  hk <- k %/% 2;
#
#  res <- diff(cumsum(c(0,x)),k)/k
#
#  #c(rep(NA, hk), res, rep(NA, hk));
#  c(rep(res[1], hk), res, rep(res[length(res)], hk));
#}

runsum <- function(x,k) {
  if ((k %% 2) == 0) { k <- k+1; }
  hk <- k %/% 2;

  res <- diff(cumsum(c(0,x)),k)

  #c(rep(NA, hk), res, rep(NA, hk));
  c(rep(res[1], hk), res, rep(res[length(res)], hk));
}

runsd <- function(x,k) {
  if ((k %% 2) == 0) { k <- k+1; }
  hk <- k %/% 2;

  a<-diff(c(0,cumsum(x^2)),k)*k
  b<-diff(c(0,cumsum(x)),k)^2
  d<-(a-b)/(k*(k-1))
  res <- sqrt(d)

  #c(rep(NA, hk), res, rep(NA, hk));
  c(rep(res[1], hk), res, rep(res[length(res)], hk));
}

runvar <- function(x,k) {
  if ((k %% 2) == 0) { k <- k+1; }
  hk <- k %/% 2;

  a<-diff(c(0,cumsum(x^2)),k)*k
  b<-diff(c(0,cumsum(x)),k)^2
  d<-(a-b)/(k*(k-1))
  res <- d

  #c(rep(NA, hk), res, rep(NA, hk));
  c(rep(res[1], hk), res, rep(res[length(res)], hk));
}

# Function to calculate a running t-statistic for samples with (unknown) unequal
# variance and size (i.e., like in Welch t-test)
runt <- function(x, k) {
  if ((k %% 2) != 0) { k <- k+1; }
  hk <- k %/% 2;

  #m <- runmean(x, k=hk)
  m <- runsum(x, k=hk)/hk
  m1 <- c(rep(NA, hk), m)
  m2 <- c(m, rep(NA, hk))
  s <- runvar(x, k=hk)
  s1 <- c(rep(NA, hk), s)
  s2 <- c(s, rep(NA, hk))

  res <- (m1 - m2) / sqrt(s1/hk + s2/hk)
  hhk <- hk %/% 2;
  res[(hhk+1):(length(res)-(hhk+1))]
}

## Originally by Elzo de Wit (does not work properly)
#runcor <- function(a,b,n) {
#  ab <- a * b
#  prod <- diff(c(0, cumsum(ab)), n)
#  #sub <- n * runmean(a,n) * runmean(b,n)
#  sub <- n * (runsum(a, k=n)/n) * (runsum(b, k=n)/n)
#  bottom <- (n-1) * runsd(a,n) * runsd(b,n)
#  (prod-sub)/bottom
#}


odd <- function(x) { x[x %% 2 == 1] }
even <- function(x) { x[x %% 2 == 0] }

readGFF <- function(file, purge=c(), nlines = NULL, skip = 0, verbose=FALSE) {
  if(is.null(nlines))
    nlines <- as.numeric(system(paste("cat \"", file, "\" | wc -l", sep=''), intern=TRUE))
              
  if(verbose)
    cat(paste("Reading in", nlines, "lines..."));
      
  what_list <- list(seqname='c', source='c', feature='c', start=0, end=0,
                    score=0.1, strand='0', frame='c', attribute='c');
        
  if (length(purge) > 0) { what_list[purge] <- list(NULL); }
          
  data <- scan(file, nmax=nlines, skip=skip, what=what_list, sep='\t', quiet=TRUE);
            
  if (length(purge) > 0) { data[purge] <- NULL }
              
  class(data) <- "data.frame";
  row.names(data) <- 1:nlines;
                  
  if(verbose)
    cat("done.\n");
                    
  data$score <- as.numeric(data$score);
  if (all(is.na(data$score))) { data$score <- 0; }

  invisible(data);
}


writeGFF <- function(data, file) {
  n <- nrow(data);
  newdata <- data.frame( seqname= rep(".", n), source = rep(".", n), feature   = rep(".", n),
                         start  = rep(0, n),   end    = rep(0, n),   score     = rep(0, n),
                         strand = rep(".", n), frame  = rep(".", n), attribute = rep(".", n) )
  for (i in colnames(data)) {
    if (i %in% names(newdata)) {
      newdata[,i] <- data[,i];
    }
  }

  write.table(newdata, file=file, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE);
}

readBED <- function(file, purge=c(), nlines = NULL, skip = 0, verbose = FALSE) {
  if(is.null(nlines))
    nlines <- as.numeric(system(paste("cat \"", file, "\" | wc -l", sep=''), intern=TRUE))
              
  if(verbose)
    cat(paste("Reading in", nlines, "lines..."));
      
  what_list <- list(seqname='c', start=0, end=0, name='c', score='c', strand='c');
        
  if (length(purge) > 0) { what_list[purge] <- list(NULL); }
          
  data <- scan(file, nmax=nlines, skip=skip, what=what_list, sep='\t', quiet=TRUE);
            
  if (length(purge) > 0) { data[purge] <- NULL }
              
  class(data) <- "data.frame";
  row.names(data) <- 1:nlines;
                  
  if(verbose)
    cat("done.\n");
                    
  data$score <- as.numeric(data$score);
  if (all(is.na(data$score))) { data$score <- 0; }

  invisible(data);
}

readBED4 <- function(file, purge=c(), nlines = NULL, skip = 0, verbose = FALSE) {
  if(is.null(nlines))
    nlines <- as.numeric(system(paste("cat \"", file, "\" | wc -l", sep=''), intern=TRUE))
              
  if(verbose)
    cat(paste("Reading in", nlines, "lines..."));
      
  what_list <- list(seqname='c', start=0, end=0, score='c');
        
  if (length(purge) > 0) { what_list[purge] <- list(NULL); }
          
  data <- scan(file, nmax=nlines, skip=skip, what=what_list, sep='\t', quiet=TRUE);
            
  if (length(purge) > 0) { data[purge] <- NULL }
              
  class(data) <- "data.frame";
  row.names(data) <- 1:nlines;
                  
  if(verbose)
    cat("done.\n");
                    
  data$score <- as.numeric(data$score);
  if (all(is.na(data$score))) { data$score <- 0; }

  invisible(data);
}


writeBED <- function(data, file, ...) {
  n <- nrow(data);
  newdata <- data.frame( seqname= rep(".", n), start = rep(0, n), end    = rep(0, n),   
                         name  = rep(".", n),  score = rep(0, n), strand = rep(".", n) )
  for (i in colnames(data)) {
    if (i %in% names(newdata)) {
      newdata[,i] <- data[,i];
    }
  }

  write.table(newdata, file=file, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE, ...);
}

writeWIG <- function(data, file) {
  n <- nrow(data);
  chroms <- unique(data$seqname);
  span <- median(data$end-data$start);
  append <- FALSE;
  for (chr in chroms) {
    write(paste("variableStep chrom=", chr, " span=", span, sep=""), file=file, append=append)
    write.table(data[data$seqname==chr, c("start", "score")], file=file, 
                quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE, append=TRUE);
    append <- TRUE;
  }
}

get_overlap <- function(regions, feats, modus=1, minOverlap=0.5) {
  region_feats = NULL;
  # We loop over regions, because we generally have
  # less regions than features.
  for (i in 1:nrow(regions)) {
    if (i %% 100 == 0) { cat(paste(i, "\r")); }

    chr_hits <- feats$seqname == regions$seqname[i];
    cont_hits <- 1; # By default, this is satisfied. Only used for modus 3.
    # Options: 1 => fully contained, 2 => partly contained,
    #          3 => fraction contained at least 'minOverlap'
    if (modus == 1) {
      start_hits <- feats$start > regions$start[i];
      end_hits   <- feats$end < regions$end[i];
    } else if (modus == 2 || modus == 3) {
      start_hits <- feats$end > regions$start[i];
      end_hits   <- feats$start < regions$end[i];

      if (modus == 3) {
        # Find position of feature relative to current region
        pend <- pmin(feats$end,   regions$end[i])
        pstart   <- pmax(feats$start, regions$start[i])
        overlap <- (pend - pstart) / (feats$end - feats$start)
        cont_hits <- overlap >= minOverlap
      }
    } else {
      error("Unknown modus: use either 1, 2 or 3")
    }
    hits <- chr_hits & start_hits & end_hits & cont_hits

    if (any(hits, na.rm=TRUE)) {
      hit_feats <- feats[which(hits),]
      feats <- feats[-which(hits),];
      hit_feats$frame = paste("Region", i);
      region_feats <- rbind(region_feats, hit_feats);
    }
  }
  cat("\n");

  list(fg=region_feats, bg=feats)
}


plot_regional_expression <- function(regions, file, delta=FALSE) {
  genes <- readGFF("/data/people/meuleman/projects/cons/LADs_fourth_strike/data/Ensembl50_Mmus_NCBI37_genes.gff")
  exons <- readGFF("/data/people/meuleman/projects/cons/LADs_fourth_strike/data/Ensembl50_Mmus_NCBI37_exons.gff")
  data_path <- "/data/datasets/Wouter/Efroni2008_AffyTiling/mm16tpmap/"; # Trailing slash!

  if (delta) {
    pdf(file=file, width=10, version="1.4")
    par(mar=c(4,4,2,2))
    par(mfrow=c(3,1))
  } else {
    pdf(file=file, width=14)
    par(mar=c(4,4,2,2))
    par(mfrow=c(2,1))
  }

  for (region in regions) {
    print(region);
    items <- unlist(strsplit(region, "[:-]"));
    chr <- items[1]; 
    start <- as.numeric(items[2]); 
    end <- as.numeric(items[3]);

    exons_chr <- exons[exons$seqname==chr,];
    genes_chr <- genes[genes$seqname==chr,];

    ES_name <- load(paste(data_path, chr, "_ES.RData", sep=""))
    ES <- get(ES_name)
    NP_name <- load(paste(data_path, chr, "_NP.RData", sep=""))
    NP <- get(NP_name)

    #ES$score[ES$score <= 0] <- NA
    #NP$score[NP$score <= 0] <- NA

    ES$score <- ES$score - 4
    NP$score <- NP$score - 4

    sel <- (NP$start > start & NP$end < end)
    ymax <- max(c(ES$score[sel], NP$score[sel]), na.rm=T);

    plot(0, type="n", xlim=c(start,end), ylim=c(-ymax/5,ymax), xlab="Genomic location", 
         ylab="signal", col=2, main=paste("ESC (", region, ")", sep=""))
    rect(ES$start, 0, ES$end, ES$score, border="NA", col=2)
    # Genes
    rect(genes_chr$start[genes_chr$strand=="+"], -1, genes_chr$end[genes_chr$strand=="+"],  -ymax/10, border=TRUE)
    rect(genes_chr$start[genes_chr$strand=="-"], -ymax/10, genes_chr$end[genes_chr$strand=="-"], -ymax/5, border=TRUE)
    # Exons
    rect(exons_chr$start[exons_chr$strand=="+"], -1, exons_chr$end[exons_chr$strand=="+"],  -ymax/10, 
      col=as.factor(exons_chr$attribute[exons_chr$strand=="+"]))
    rect(exons_chr$start[exons_chr$strand=="-"], -ymax/10, exons_chr$end[exons_chr$strand=="-"], -ymax/5, 
      col=as.factor(exons_chr$attribute[exons_chr$strand=="-"]))

    plot(0, type="n", xlim=c(start,end), ylim=c(-ymax/5,ymax), xlab="Genomic location", 
         ylab="signal", col=3, main=paste("NPC (", region, ")", sep=""))
    rect(NP$start, 0, NP$end, NP$score, border="NA", col=3)
    # Genes
    rect(genes_chr$start[genes_chr$strand=="+"], -1, genes_chr$end[genes_chr$strand=="+"],  -ymax/10, border=TRUE)
    rect(genes_chr$start[genes_chr$strand=="-"], -ymax/10, genes_chr$end[genes_chr$strand=="-"], -ymax/5, border=TRUE)
    # Exons
    rect(exons_chr$start[exons_chr$strand=="+"], -1, exons_chr$end[exons_chr$strand=="+"],  -ymax/10, 
      col=as.factor(exons_chr$attribute[exons_chr$strand=="+"]))
    rect(exons_chr$start[exons_chr$strand=="-"], -ymax/10, exons_chr$end[exons_chr$strand=="-"], -ymax/5, 
      col=as.factor(exons_chr$attribute[exons_chr$strand=="-"]))

    if (delta) {
      ratio_name <- load(paste(data_path, chr, ".RData", sep=""))
      ratio <- get(ratio_name)
      dScore <- ratio$score

      nonNA  <- (!is.na(ES$score) & !is.na(NP$score) & !is.na(dScore));
      ES <- ES[nonNA,]
      NP <- NP[nonNA,]
      dScore <- dScore[nonNA]

      # Used for running medians
      ES$start <- NP$start <- runmed(NP$start, k=49)
      ES$end <- NP$end <- runmed(NP$end, k=49)
      ES$score <- runmed(ES$score, k=49)
      NP$score <- runmed(NP$score, k=49)
      dScore <- runmed(dScore, k=49);

      #dScore <- log2(ES$score/NP$score);
      sel <- (NP$start > start & NP$end < end)

      yrange <- range(c(0, dScore[sel]), finite=TRUE);
      plot(0, type="n", xlim=c(start,end), ylim=yrange, xlab="Genomic location", 
           ylab="signal", main=paste("Z-score (", region, ")", sep=""))

      use_alpha <- FALSE;
      if (use_alpha) {
        signal <- ES$score+NP$score
        alphas <- (signal - min(signal)) / (max(signal) - min(signal))
        rect(ES$start[dScore>0], 0, ES$end[dScore>0], dScore[dScore>0], 
             border="NA", col=rgb(1,0,0,alphas[dScore>0]))
        rect(NP$start[dScore<0], dScore[dScore<0], NP$end[dScore<0], 0, 
             border="NA", col=rgb(0,1,0,alphas[dScore<0]))
      } else {
        rect(ES$start[dScore>0], 0, ES$end[dScore>0], dScore[dScore>0], 
             border="NA", col=rgb(1,0,0))
        rect(NP$start[dScore<0], dScore[dScore<0], NP$end[dScore<0], 0, 
             border="NA", col=rgb(0,1,0))
      }
    }
  }
  
  dev.off()
}


CGintvlOverlap_Wouter <- function (domains, features, remapIndices = TRUE, fun=mean, ...) {
    allseq = unique(c(domains$seqname, features$seqname))
    res = vector(mode = "list", length = length(allseq))
    names(res) = allseq
    warnings <- c()
    for (chr in allseq) {
        print(chr)
        d <- domains[domains$seqname == chr, ]
        fOrder <- which(features$seqname == chr)
        f <- features[fOrder, ]
        res[[chr]] <- .domOverlap(d$start, d$end, f$start, f$end, 
            ...)
        warnings <- c(warnings, res[[chr]]$warn)
        res[[chr]] <- res[[chr]]$overlap

        #res[[chr]] <- lapply(res[[chr]], function(v) median(f$score[v], na.rm=T))
        res[[chr]] <- lapply(res[[chr]], function(v) fun(f$score[v]))
        #res[[chr]] <- lapply(res[[chr]], function(v) f$score[v])
        #res[[chr]] <- lapply(res[[chr]], function(v) min(f$score[v]))
#        if (remapIndices) 
#            res[[chr]] <- lapply(res[[chr]], function(v) fOrder[v])
    }
    res
}

# For every (extended) domain, the function finds the features contained in it.
# Then, the index of the nearest neighbour feature (based on domain and feature
# centers) is returned. If no features are found, NA is returned.
CGintvlOverlap_Wouter_closest <- function (domains, features, ...) {
    allseq = unique(c(domains$seqname, features$seqname))
    res = vector(mode = "list", length = length(allseq))
    names(res) = allseq
    warnings <- c()
    for (chr in allseq) {
        print(chr)
        d <- domains[domains$seqname == chr, ]
        fOrder <- which(features$seqname == chr)
        f <- features[fOrder, ]
        dctr <- (d$start+d$end)/2
        fctr <- (f$start+f$end)/2
        res[[chr]] <- .domOverlap(d$start, d$end, f$start, f$end, ...)
        warnings <- c(warnings, res[[chr]]$warn)
        res[[chr]] <- res[[chr]]$overlap

        di <- 0;
        res[[chr]] <- lapply(res[[chr]], function(v) {
            di <<- di + 1;
            if (length(v) > 0)
              fOrder[v[which.min(abs(fctr[v] - dctr[di]))]]
            else
              NA
        })
    }
    res
}


fix_overlap <- function(x.in) {
  # create matrix to determine queue size (overlap)
  x.q <- rbind(cbind(x.in$start, 1, 1:nrow(x.in)), 
               cbind(x.in$end+1,  -1, 1:nrow(x.in)))
  # sort
  x.q <- x.q[order(x.q[,1], x.q[,2]),]
  x.q <- cbind(x.q, queue=cumsum(x.q[,2]))

  end   <- which(x.q[,4] == 0);
  start <- c(1, end+1);
  start <- start[-length(start)];

  # This removes all overlapping features!
  #diff_res  <- apply(cbind(start, end), 1, diff)
  #diff_hits <- which(diff_res > 1)
  #diff_sel  <- unlist(apply(cbind(start[diff_hits], 
  #                          end[diff_hits]), 
  #                    1, function(x) x[1]:x[2]))
  #x.in[-unique(x.q[diff_sel,3]),]

  data.frame(start=x.q[start,1], end=x.q[end,1]-1)
}

plot_delta_dom_sizes <- function(delta_doms_pos, delta_doms_neg, chr_sizes, file="delta_doms_sizes.pdf") {
  pdf(file=file)

  delta_doms  <- rbind(delta_doms_pos, delta_doms_neg);
  sizes       <- delta_doms$end-delta_doms$start+1

  empty_sizes <- chr_sizes; empty_sizes[names(chr_sizes)] <- 0;

  delta_sizes <- pos_delta_sizes <- neg_delta_sizes <- empty_sizes
  tmp_sizes <- tapply(sizes, delta_doms$seqname, function(x) sum(x))
  delta_sizes[names(tmp_sizes)] <- tmp_sizes

  print(cbind(median(sizes), length(sizes)));
  plot(dens <- density(log10(sizes)), lwd=2,
       main="Size of delta Lamin domains (+ & -)", xlab="log10(bp)", ylab="Density")
  abline(v=median(log10(sizes)), lwd=2, col=2)
  
  barplot(round(delta_sizes/chr_sizes*100, 2), las=2, ylab="% of chromosome size", 
          main=paste("delta Lamin fractions (+ & -)\nGenome-wide: ", 
               round(sum(as.numeric(delta_sizes))/sum(as.numeric(chr_sizes))*100, 2), "%", sep=""))
  
  delta_sizes <- chr_sizes
  tmp_sizes <- tapply(sizes, delta_doms$seqname, function(x) sum(x))
  delta_sizes[names(delta_sizes)] <- 0;
  delta_sizes[names(tmp_sizes)] <- tmp_sizes

  pos_sizes       <- delta_doms_pos$end-delta_doms_pos$start+1
  pos_delta_sizes <- empty_sizes
  tmp_sizes <- tapply(pos_sizes, delta_doms_pos$seqname, function(x) sum(x))
  pos_delta_sizes[names(tmp_sizes)] <- tmp_sizes
  print(cbind(median(pos_sizes), length(pos_sizes)));
  plot(dens <- density(log10(pos_sizes)), lwd=2,
       main="Size of delta Lamin domains (+)", xlab="log10(bp)", ylab="Density")
  abline(v=median(log10(pos_sizes)), lwd=2, col=2)
  
  pos_delta_frac <- round(pos_delta_sizes/chr_sizes*100, 2)
  barplot(pos_delta_frac, las=2, ylab="% of chromosome size", col="blue",
          main=paste("delta Lamin fractions (+)\nGenome-wide: ", 
               round(sum(as.numeric(pos_delta_sizes))/sum(as.numeric(chr_sizes))*100, 2), "%", sep=""))
  
  neg_sizes       <- delta_doms_neg$end-delta_doms_neg$start+1
  neg_delta_sizes <- empty_sizes
  tmp_sizes <- tapply(neg_sizes, delta_doms_neg$seqname, function(x) sum(x))
  neg_delta_sizes[names(tmp_sizes)] <- tmp_sizes
  print(cbind(median(neg_sizes), length(neg_sizes)));
  plot(dens <- density(log10(neg_sizes)), lwd=2,
       main="Size of delta Lamin domains (-)", xlab="log10(bp)", ylab="Density")
  abline(v=median(log10(neg_sizes)), lwd=2, col=2)
  
  neg_delta_frac <- round(neg_delta_sizes/chr_sizes*100, 2)
  barplot(neg_delta_frac, las=2, ylab="% of chromosome size", col="orange",
          main=paste("delta Lamin fractions (-)\nGenome-wide: ", 
               round(sum(as.numeric(neg_delta_sizes))/sum(as.numeric(chr_sizes))*100, 2), "%", sep=""))
  
  pos_neg_mat <- rbind(pos_delta_frac, neg_delta_frac)
  barplot(pos_neg_mat, las=2, ylab="% of chromosome size", col=c("blue","orange"),
          main=paste("delta Lamin fractions (+ & -)\nGenome-wide: ", 
               round(sum(as.numeric(delta_sizes))/sum(as.numeric(chr_sizes))*100, 2), "%", sep=""))
  legend("topright", "(x,y)", fill=c("blue","orange"), legend=c("+", "-"))
  
  dev.off()
}

plot_detail <- function(ES, NP, NP_ES, delta_doms_pos, delta_doms_neg, chr, from, to, id="default", path="plots", pdf=FALSE, lwd=1) {
  if (pdf) {
    pdf(file=paste(path, "/detail_", chr, "_", from, "-", to, "_", id, ".pdf", sep=""), width=15)
  } else {
#    png(file=paste(path, "/detail_", chr, "_", from, "-", to, "_", id, ".png", sep=""), width=960)
  }

  ES_chr <- ES[ES$seqname==chr,]
  NP_chr <- NP[ES$seqname==chr,]
  NP_ES_chr <- NP_ES[ES$seqname==chr,]
  delta_doms_neg_chr <- delta_doms_neg[delta_doms_neg$seqname==chr,]
  delta_doms_pos_chr <- delta_doms_pos[delta_doms_pos$seqname==chr,]
  
  ext <- (to - from)
  print(paste("Plotting region ", chr, ":", from-ext, "-", to+ext, sep=""))

  par(mfrow=c(4,1), mar=c(0,2,0,1))
  plot(ES_chr$start, runmed(ES_chr$score, k=5), type="h", col="orange", xlim=c(from-ext, to+ext), xaxt="n", xaxs="i", lwd=lwd)
  abline(v=c(from,to), lty=2)
  plot(NP_chr$start, runmed(NP_chr$score, k=5), type="h", col="blue", xlim=c(from-ext, to+ext), xaxt="n", xaxs="i", lwd=lwd)
  abline(v=c(from,to), lty=2)
  plot(NP_ES_chr$start, runmed(NP_ES_chr$score, k=5), type="h", col=ifelse(runmed(NP_ES_chr$score, k=5) > 0, "blue", "orange"), xlim=c(from-ext, to+ext), xaxt="n", xaxs="i", lwd=lwd)
  abline(v=c(from,to), lty=2)
  plot(0, type="n", xlim=c(from-ext, to+ext), ylim=c(0, 2), yaxt="n", axes=FALSE, xlab="", ylab="", xaxs="i")
  par(mar=c(8,2,0,1), new=TRUE)
  plot(0, type="n", xlim=c(from-ext, to+ext), ylim=c(0, 2), yaxt="n", xaxs="i", xlab=paste("Chromosome", sub("chr", "", chr), "(bp)"))
  if (!is.null(delta_doms_neg_chr))
    rect(delta_doms_neg_chr$start, 1.5, delta_doms_neg_chr$end, 1.7, col="grey")
  if (!is.null(delta_doms_pos_chr))
    rect(delta_doms_pos_chr$start, 1.8, delta_doms_pos_chr$end, 2, col="grey")

#  dev.off()
}


plot_three_singles <- function(ES, NP, EF, ES_doms, NP_doms, EF_doms, id="default", path="plots", pdf=FALSE) {
  chroms <- unique(ES$seqname)
  for (i in chroms) {
    print(i)

    if (pdf) {
      pdf(file=paste(path, "/three_singles_", id, "_", i, ".pdf", sep=""), width=15)
    } else {
      png(file=paste(path, "/three_singles_", id, "_", i, ".png", sep=""), width=960)
    }

    par(mfrow=c(3,1))

    ylim <- quantile(ES$score[ES$seqname == i], probs=c(0.01,0.99))
    plot(ES$start[ES$seqname == i], runmed(ES$score[ES$seqname == i], 49), type="h", bty="n", 
         xlab="genomic location (bp)", ylab="LaminB1 interaction", 
         ylim=ylim, main="ES cells", col="orange")
    if (!is.null(ES_doms) && (nrow(ES_doms[ES_doms$seqname == i,]) > 0)) 
        rect(ES_doms$start[ES_doms$seqname == i], ylim[1], ES_doms$end[ES_doms$seqname == i], ylim[1]+0.2, col = "grey");

    ylim <- quantile(NP$score[NP$seqname == i], probs=c(0.01,0.99))
    plot(NP$start[NP$seqname == i], runmed(NP$score[NP$seqname == i], 49), type="h", bty="n", 
         xlab="genomic location (bp)", ylab="LaminB1 interaction", 
         ylim=ylim, main="NP cells", col="blue")
    if (!is.null(NP_doms) && (nrow(NP_doms[NP_doms$seqname == i,]) > 0)) 
      rect(NP_doms$start[NP_doms$seqname == i], ylim[1], NP_doms$end[NP_doms$seqname == i], ylim[1]+0.2, col = "grey");

    ylim <- quantile(EF$score[EF$seqname == i], probs=c(0.01,0.99))
    plot(EF$start[EF$seqname == i], runmed(EF$score[EF$seqname == i], 49), type="h", bty="n", 
         xlab="genomic location (bp)", ylab="LaminB1 interaction", 
         ylim=ylim, main="EF cells", col="magenta")
    if (!is.null(EF_doms) && (nrow(EF_doms[EF_doms$seqname == i,]) > 0)) 
      rect(EF_doms$start[EF_doms$seqname == i], ylim[1], EF_doms$end[EF_doms$seqname == i], ylim[1]+0.2, col = "grey");

    dev.off()
  }
}

plot_detail_specific <- function(signal1, signal2, delta, Entrez, sel) {
  chr <- min(sel)
  plot_detail(signal1, signal2, delta, NULL, NULL, chr=Entrez$seqname[chr], 
    from=Entrez$start[min(sel)], to=Entrez$end[max(sel)])
  rect(Entrez$start[Entrez$seqname==Entrez$seqname[chr] & Entrez$strand=="+"], 1,
       Entrez$end[Entrez$seqname==Entrez$seqname[chr] & Entrez$strand=="+"], 1.2, col=2)
  rect(Entrez$start[Entrez$seqname==Entrez$seqname[chr] & Entrez$strand=="-"], 0.7,
       Entrez$end[Entrez$seqname==Entrez$seqname[chr] & Entrez$strand=="-"], 0.9, col=2)
}

plot_detail_specific_Ensembl <- function(signal1, signal2, delta, Entrez, Ensembl, sel) {
  chr <- min(sel)
  plot_detail(signal1, signal2, delta, NULL, NULL, chr=Entrez$seqname[chr], 
    from=Entrez$start[min(sel)], to=Entrez$end[max(sel)])

  rect(Entrez$start[Entrez$seqname==Entrez$seqname[chr] & Entrez$strand=="+"], 1,
       Entrez$end[Entrez$seqname==Entrez$seqname[chr] & Entrez$strand=="+"], 1.2, col=2)
  rect(Entrez$start[Entrez$seqname==Entrez$seqname[chr] & Entrez$strand=="-"], 0.7,
       Entrez$end[Entrez$seqname==Entrez$seqname[chr] & Entrez$strand=="-"], 0.9, col=2)

  rect(Ensembl$start[Ensembl$seqname==Entrez$seqname[chr] & Ensembl$strand=="+"], 0.4,
       Ensembl$end[Ensembl$seqname==Entrez$seqname[chr] & Ensembl$strand=="+"], 0.6, col=rgb(0,1,0,0.3))
  rect(Ensembl$start[Ensembl$seqname==Entrez$seqname[chr] & Ensembl$strand=="-"], 0.3,
       Ensembl$end[Ensembl$seqname==Entrez$seqname[chr] & Ensembl$strand=="-"], 0.5, col=rgb(0,1,0,0.3))
}

plot_single_and_comb <- function(ES, NP, NP_ES, ES_doms, NP_doms, delta_doms_pos, delta_doms_neg, comb="delta", id="default", path="plots", pdf=FALSE) {
  chroms <- unique(ES$seqname)
  for (i in chroms) {
    print(i)

    if (pdf) {
      pdf(file=paste(path, "/single_and_", comb, "_", id, "_", i, ".pdf", sep=""), width=15)
    } else {
      png(file=paste(path, "/single_and_", comb, "_", id, "_", i, ".png", sep=""), width=960)
    }

    par(mfrow=c(3,1))

    ylim <- quantile(ES$score[ES$seqname == i], probs=c(0.01,0.99))
    plot(ES$start[ES$seqname == i], runmed(ES$score[ES$seqname == i], 49), type="h", bty="n", 
         xlab="genomic location (bp)", ylab="LaminB1 interaction", 
         ylim=ylim, main="ES cells", col="orange")
    if (!is.null(ES_doms) && (nrow(ES_doms[ES_doms$seqname == i,]) > 0)) 
        rect(ES_doms$start[ES_doms$seqname == i], ylim[1], ES_doms$end[ES_doms$seqname == i], ylim[1]+0.2, col = "grey");

    ylim <- quantile(NP$score[ES$seqname == i], probs=c(0.01,0.99))
    plot(NP$start[NP$seqname == i], runmed(NP$score[NP$seqname == i], 49), type="h", bty="n", 
         xlab="genomic location (bp)", ylab="LaminB1 interaction", 
         ylim=ylim, main="NP cells", col="blue")
    if (!is.null(NP_doms) && (nrow(NP_doms[NP_doms$seqname == i,]) > 0)) 
      rect(NP_doms$start[NP_doms$seqname == i], ylim[1], NP_doms$end[NP_doms$seqname == i], ylim[1]+0.2, col = "grey");

    ylim <- quantile(NP_ES$score[ES$seqname == i], probs=c(0.01,0.99))
    delta_signal <- runmed(NP_ES$score[NP_ES$seqname == i], 49);
    delta_col <- ifelse(delta_signal > 0, "blue", "orange")
    plot(NP_ES$start[NP_ES$seqname == i], delta_signal, type="h", bty="n", 
         xlab="genomic location (bp)", ylab="LaminB1 interaction", 
         ylim=ylim, main=paste(comb, "Lamin"), col=delta_col)

    if (nrow(delta_doms_pos[delta_doms_pos$seqname == i,]) > 0) {
      rect(delta_doms_pos$start[delta_doms_pos$seqname == i], ylim[2], 
           delta_doms_pos$end[delta_doms_pos$seqname == i], ylim[2]+0.2, col = "grey");
    }
    if (nrow(delta_doms_neg[delta_doms_neg$seqname == i,]) > 0) {
      rect(delta_doms_neg$start[delta_doms_neg$seqname == i], ylim[1], 
           delta_doms_neg$end[delta_doms_neg$seqname == i], ylim[1]+0.2, col = "grey");
    }

    dev.off()
  }
}

### OK (OLD)
#subAlign <- function(aln,idx){
#  if (is.logical(idx)) idx <- which(idx);
#
#  newIntvlIdx <- rep(NA, nrow(aln$relDom))
#  newIntvlIdx[idx] <- 1:length(idx)
#
#  alnIdx <- (aln$aln$relDomIdx %in% rownames(aln$relDom)[idx])
#  aln$aln <- aln$aln[alnIdx,]
#  aln$aln$relDomIdx <- newIntvlIdx[aln$aln$relDomIdx];
#
#  aln$relDom <- aln$relDom[idx,]
#
#  aln
#}

subAlign <- function(aln,idx){
  if (is.logical(idx)) idx <- which(idx);

  newIntvlIdx <- rep(NA, nrow(aln$relDom))
  newIntvlIdx[idx] <- 1:length(idx)

  aln$relDom <- aln$relDom[idx,]
  aln$aln <- aln$aln[aln$aln$relDomIdx %in% idx,]
  aln$aln$relDomIdx <- newIntvlIdx[aln$aln$relDomIdx];

  aln
}

plot_profile <- function(GFF, LADs=NULL, chr="chr18", col="black", main="", type="h", 
                         plot_LADs=TRUE, plot_stats=FALSE, plot_ylab=TRUE, k=49, ylim=NULL, sample=FALSE, ...) {
  if (is.null(ylim)) ylim <- quantile(GFF$score[GFF$seqname == chr], probs=c(0.01,0.99))
  pos <- GFF$start[GFF$seqname == chr];
  scores <- runmed(GFF$score[GFF$seqname == chr], k);
  if (sample) {
    sel <- sample(1:length(pos), round(length(pos)/5))
  } else {
    sel <- 1:length(pos);
  }
  if (length(col) > 1) col <- col[sel];

  plot(pos[sel], scores[sel], type=type, bty="n",
       ylab=ifelse(plot_ylab, "LaminB1 interaction (log2)", ""), xaxt="n", ylim=ylim, col=col, ...)
  mtext(main, 3, -1, cex=1.5)
  if (!is.null(LADs) && (nrow(LADs[LADs$seqname == chr,]) > 0)) {
    if (plot_LADs) {
      rect(LADs$start[LADs$seqname == chr], ylim[1], LADs$end[LADs$seqname == chr], ylim[1]+0.4, 
           col = grey(0.5), border=NA);
    }
    if (plot_stats) {
      med <- median(LADs$end-LADs$start+1)
      legend("topright", "(x,y)", bty="n", inset=c(0.1,-0.05),
             legend=c(paste(nrow(LADs), "LADs"), paste("Median size: ", med %/% 1e3, "kb", sep="")))
    }
  }
}

## Detailed profile
plot_binding_detail <- function(chr, from, to, lwd=1, extend=0.5) {
  par(mfrow=c(4,1), mar=c(0,5,0,1))
  i <- chr
 
  from <- from - extend*(to - from)
  to <- to + extend*(to - from)
 
  ylim <- quantile(ES$score[ES$seqname == i], probs=c(0.05,0.95))
    plot(ES$start[ES$seqname == i], runmed(ES$score[ES$seqname == i], 9), type="h", bty="n", lwd=lwd, lend=2,
         ylab="", xaxt="n", ylim=ylim, col="orange", xlim=c(from,to), cex.axis=1.5)
  #mtext("Embryonic stem cells", 3, -1)
  
  ylim <- quantile(NP$score[NP$seqname == i], probs=c(0.05,0.95))
  plot(NP$start[NP$seqname == i], runmed(NP$score[NP$seqname == i], 9), type="h", bty="n", lwd=lwd, lend=2,
       ylab="", xaxt="n", ylim=ylim, col="blue", xlim=c(from,to), cex.axis=1.5)
  #mtext("Neuro-precursor cells", 3, -1)
  
  ylim <- quantile(AC$score[AC$seqname == i], probs=c(0.05,0.95))
  plot(AC$start[AC$seqname == i], runmed(AC$score[AC$seqname == i], 9), type="h", bty="n", lwd=lwd, lend=2,
       ylab="", xaxt="n", ylim=ylim, col="magenta", xlim=c(from,to), cex.axis=1.5)
  #mtext("Astrocyte cells", 3, -1)
  
  mtext("LaminB1 interaction (log2)", 2, -2, cex=1.5, at=0.62, outer=TRUE)

  par(mar=c(5,5,0,1))
  plot(AC$start[AC$seqname == i], runmed(AC$score[AC$seqname == i], 49), 
       type="n", bty="n", ylim=c(0.4, 1), xlab="", ylab="", xlim=c(from,to), yaxt="n", cex.axis=1.5, cex.lab=1.5)
  mtext(paste("Genomic location (bp, ", chr, ")", sep=""), 1, 3, cex=1.5)

  rect(Entrez_full$start[Entrez_full$seqname == i & Entrez_full$strand == "+"], 0.77, Entrez_full$end[Entrez_full$seqname == i & Entrez_full$strand == "+"], 0.78, col=1);
  rect(Entrez_full$start[Entrez_full$seqname == i & Entrez_full$strand == "-"], 0.57, Entrez_full$end[Entrez_full$seqname == i & Entrez_full$strand == "-"], 0.58, col=1);

  text(Entrez_full$start[Entrez_full$seqname == i & Entrez_full$strand == "+"], 0.87, labels=Entrez_full$frame[Entrez_full$seqname == i & Entrez_full$strand == "+"], adj=c(0,0), col=1, cex=1.5)
  text(Entrez_full$end[Entrez_full$seqname == i & Entrez_full$strand == "-"], 0.67, labels=Entrez_full$frame[Entrez_full$seqname == i & Entrez_full$strand == "-"], adj=c(1,0), col=1, cex=1.5)
  
  rect(UTRs_full$start[UTRs_full$seqname == i & UTRs_full$strand == "+"], 0.74, UTRs_full$end[UTRs_full$seqname == i & UTRs_full$strand == "+"], 0.81, col=1);
  rect(UTRs_full$start[UTRs_full$seqname == i & UTRs_full$strand == "-"], 0.54, UTRs_full$end[UTRs_full$seqname == i & UTRs_full$strand == "-"], 0.61, col=1);
  
  rect(CDSs_full$start[CDSs_full$seqname == i & CDSs_full$strand == "+"], 0.7, CDSs_full$end[CDSs_full$seqname == i & CDSs_full$strand == "+"], 0.85, col=1);
  rect(CDSs_full$start[CDSs_full$seqname == i & CDSs_full$strand == "-"], 0.5, CDSs_full$end[CDSs_full$seqname == i & CDSs_full$strand == "-"], 0.65, col=1);
}

## Detailed profile, including H3K36me3 data
plot_binding_detail_plus_K36 <- function(chr, from, to, lwd=1, extend=0.5, ylim_LB1=NULL, ylim_K36=NULL, 
                                         axes_LB1=TRUE, axes_K36=TRUE, emph_used_genes=FALSE, cols=c("orange", "blue", "magenta")) {
  par(mfrow=c(3,1), oma=c(5,0,0,0))
  i <- chr
 
  load(paste("../../data/external_data/Mikkelsen_Nat07/H3K36me3_data/mm9/", chr, ".ES.K36.RData", sep=""))
  names(K36) <- c("seqname", "start", "end", "score")
  #K36$score <- log2(K36$score)
  K36 <- K36[K36$score>2,]
  ES_K36 <- K36;

  load(paste("../../data/external_data/Mikkelsen_Nat07/H3K36me3_data/mm9/", chr, ".NP.K36.RData", sep=""))
  names(K36) <- c("seqname", "start", "end", "score")
  #K36$score <- log2(K36$score)
  K36 <- K36[K36$score>2,]
  NP_K36 <- K36;

  if (is.null(ylim_LB1)) 
    ylim_LB1 <- quantile(c(ES$score[ES$seqname == i], NP$score[NP$seqname == i], AC$score[AC$seqname == i]), probs=c(0.05,0.95))
  if (is.null(ylim_K36)) 
    ylim_K36 <- quantile(c(ES_K36$score[ES_K36$seqname == i], NP_K36$score[NP_K36$seqname == i]), probs=c(0.01,0.99))

  from <- from - extend*(to - from)
  to <- to + extend*(to - from)
 
  if (axes_LB1) {
    #yaxt_LB1="s"; ylab_LB1=expression(log[2](LaminB1~interaction));
    yaxt_LB1="s"; ylab_LB1=""; # We want to plot one large label over all 3 profiles...
  } else {
    yaxt_LB1="n"; ylab_LB1="";
  }

  par(mar=c(5,5,3,5))
  plot(ES$start[ES$seqname == i], runmed(ES$score[ES$seqname == i], 9), type="h", bty="n", lwd=lwd, lend=2,
       xlab="", xaxt="n", ylim=ylim_LB1, col=cols[1], xlim=c(from,to), cex.axis=1.5, yaxt=yaxt_LB1, ylab=ylab_LB1, cex.lab=2, cex.axis=2)
  #mtext("Embryonic stem cells", 3, -1)
  par(mar=c(0,5,15,5), new=TRUE)
  plot(0, type="n", bty="n", xlab="", ylab="", xaxt="n", yaxt="n", xlim=c(from, to), ylim=ylim_K36, yaxs="i")
  #rect(0, 0, 1e9, 100, col="orange", border=NA)
  lines(ES_K36$start[ES_K36$seqname == i], ES_K36$score[ES_K36$seqname == i], type="h", bty="n", lend=2, col="black")
  lines(c(ylim_K36[1], 1e9), rep(ylim_K36[1], 2), type="l", col=cols[1], lwd=3);
  if (axes_K36) {
    axis(4, cex.axis=2)
    mtext("H3K36me3", 4, line=3, cex=1.5)
  }
  
  par(mar=c(5,5,3,5))
  plot(NP$start[NP$seqname == i], runmed(NP$score[NP$seqname == i], 9), type="h", bty="n", lwd=lwd, lend=2,
       xlab="", xaxt="n", ylim=ylim_LB1, col=cols[2], xlim=c(from,to), cex.axis=1.5, yaxt=yaxt_LB1, ylab=ylab_LB1, cex.lab=2, cex.axis=2)
  #mtext("Neuro-precursor cells", 3, -1)
  par(mar=c(0,5,15,5), new=TRUE)
  plot(0, type="n", bty="n", xlab="", ylab="", xaxt="n", yaxt="n", xlim=c(from, to), ylim=ylim_K36, yaxs="i")
  #rect(0, 0, 1e9, 100, col="blue", border=NA)
  lines(NP_K36$start[NP_K36$seqname == i], NP_K36$score[NP_K36$seqname == i], type="h", bty="n", lend=2, col="black")
  lines(c(ylim_K36[1], 1e9), rep(ylim_K36[1], 2), type="l", col=cols[2], lwd=3);
  if (axes_K36) {
    axis(4, cex.axis=2)
    mtext("H3K36me3", 4, line=3, cex=1.5)
  }
  
  par(mar=c(5,5,3,5))
  plot(AC$start[AC$seqname == i], runmed(AC$score[AC$seqname == i], 9), type="h", bty="n", lwd=lwd, lend=2,
       xlab="", xaxt="n", ylim=ylim_LB1, col=cols[3], xlim=c(from,to), cex.axis=1.5, yaxt=yaxt_LB1, ylab=ylab_LB1, cex.lab=2, cex.axis=2)
  #mtext("Astrocyte cells", 3, -1)
  
  #yaxt_LB1="s"; ylab_LB1=expression(log[2](LaminB1~interaction));

  mtext(expression(log[2](LaminB1~interaction)), 2, -2, cex=1.5, at=0.55, outer=TRUE)
  #mtext("H3K36me3 (ChIP-seq data)", 4, -2, cex=1.5, at=0.53, outer=TRUE)

  #par(mar=c(5,5,0,5))
  par(new=TRUE)
  par(mar=c(0,5,3,5))
  options("scipen"=5)
  plot(AC$start[AC$seqname == i], runmed(AC$score[AC$seqname == i], 49), 
       type="n", bty="n", ylim=c(0.4, 2.1), xlab="", ylab="", xlim=c(from,to), yaxt="n", cex.axis=2)
  options("scipen"=0)
  mtext(paste("chromosomal location (bp, ", chr, ")", sep=""), 1, 3.5, cex=1.5)

  rect(Entrez_full$start[Entrez_full$seqname == i & Entrez_full$strand == "+"], 0.77, Entrez_full$end[Entrez_full$seqname == i & Entrez_full$strand == "+"], 0.78, border=NA, col=1);
  rect(Entrez_full$start[Entrez_full$seqname == i & Entrez_full$strand == "-"], 0.57, Entrez_full$end[Entrez_full$seqname == i & Entrez_full$strand == "-"], 0.58, border=NA, col=1);

  text(Entrez_full$start[Entrez_full$seqname == i & Entrez_full$strand == "+"], 0.88, labels=Entrez_full$frame[Entrez_full$seqname == i & Entrez_full$strand == "+"], adj=c(0,0), col=1, cex=2)
  text(Entrez_full$end[Entrez_full$seqname == i & Entrez_full$strand == "-"], 0.68, labels=Entrez_full$frame[Entrez_full$seqname == i & Entrez_full$strand == "-"], adj=c(1,0), col=1, cex=2)
  
  rect(UTRs_full$start[UTRs_full$seqname == i & UTRs_full$strand == "+"], 0.74, UTRs_full$end[UTRs_full$seqname == i & UTRs_full$strand == "+"], 0.81, border=NA, col=1);
  rect(UTRs_full$start[UTRs_full$seqname == i & UTRs_full$strand == "-"], 0.54, UTRs_full$end[UTRs_full$seqname == i & UTRs_full$strand == "-"], 0.61, border=NA, col=1);
  
  rect(CDSs_full$start[CDSs_full$seqname == i & CDSs_full$strand == "+"], 0.7, CDSs_full$end[CDSs_full$seqname == i & CDSs_full$strand == "+"], 0.85, border=NA, col=1);
  rect(CDSs_full$start[CDSs_full$seqname == i & CDSs_full$strand == "-"], 0.5, CDSs_full$end[CDSs_full$seqname == i & CDSs_full$strand == "-"], 0.65, border=NA, col=1);

  # Again, only genes used in analysis, in red.
  if (emph_used_genes) {
    rect(Entrez$start[Entrez$seqname == i & Entrez$strand == "+"], 0.77, Entrez$end[Entrez$seqname == i & Entrez$strand == "+"], 0.78, border=NA, col=2);
    rect(Entrez$start[Entrez$seqname == i & Entrez$strand == "-"], 0.57, Entrez$end[Entrez$seqname == i & Entrez$strand == "-"], 0.58, border=NA, col=2);

    text(Entrez$start[Entrez$seqname == i & Entrez$strand == "+"], 0.88, labels=Entrez$frame[Entrez$seqname == i & Entrez$strand == "+"], adj=c(0,0), col=2, cex=2)
    text(Entrez$end[Entrez$seqname == i & Entrez$strand == "-"], 0.68, labels=Entrez$frame[Entrez$seqname == i & Entrez$strand == "-"], adj=c(1,0), col=2, cex=2)
  
    rect(UTRs$start[UTRs$seqname == i & UTRs$strand == "+"], 0.74, UTRs$end[UTRs$seqname == i & UTRs$strand == "+"], 0.81, border=NA, col=2);
    rect(UTRs$start[UTRs$seqname == i & UTRs$strand == "-"], 0.54, UTRs$end[UTRs$seqname == i & UTRs$strand == "-"], 0.61, border=NA, col=2);
  
    rect(CDSs$start[CDSs$seqname == i & CDSs$strand == "+"], 0.7, CDSs$end[CDSs$seqname == i & CDSs$strand == "+"], 0.85, border=NA, col=2);
    rect(CDSs$start[CDSs$seqname == i & CDSs$strand == "-"], 0.5, CDSs$end[CDSs$seqname == i & CDSs$strand == "-"], 0.65, border=NA, col=2);
  }

  #sel_genes <- Entrez_full$frame[Entrez_full$seqname == i & Entrez_full$start > from & Entrez_full$end < to]
  sel_genes <- Entrez$frame[Entrez$seqname == i & Entrez$start > from & Entrez$end < to]
  invisible(sel_genes)
}

## Detailed profile, including H3K36me3 data
plot_binding_detail_plus_K36_NP_AC <- function(signal1=NP, signal2=AC, chr, from, to, lwd=1, extend=0.5, ylim_LB1=NULL, ylim_K36=NULL, 
                                         axes_LB1=TRUE, axes_K36=TRUE, emph_used_genes=FALSE, cols=c("orange", "blue", "magenta"),
                                         ylab_text=expression(log[2](LaminB1~interaction))) {
  par(mfrow=c(2,1), oma=c(5,0,0,0))
  i <- chr
 
  load(paste("../../data/external_data/Mikkelsen_Nat07/H3K36me3_data/mm9/", chr, ".NP.K36.RData", sep=""))
  names(K36) <- c("seqname", "start", "end", "score")
  #K36$score <- log2(K36$score)
  K36 <- K36[K36$score>2,]
  signal1_K36 <- K36;

  if (is.null(ylim_LB1)) 
    ylim_LB1 <- quantile(c(signal1$score[signal1$seqname == i], signal2$score[signal2$seqname == i]), probs=c(0.05,0.95))
  if (is.null(ylim_K36)) 
    ylim_K36 <- quantile(c(signal1_K36$score[signal1_K36$seqname == i]), probs=c(0.01,0.99))

  from <- from - extend*(to - from)
  to <- to + extend*(to - from)
 
  if (axes_LB1) {
    #yaxt_LB1="s"; ylab_LB1=expression(log[2](LaminB1~interaction));
    yaxt_LB1="s"; ylab_LB1=""; # We want to plot one large label over all 3 profiles...
  } else {
    yaxt_LB1="n"; ylab_LB1="";
  }

  par(mar=c(5,5,3,5))
  plot(signal1$start[signal1$seqname == i], runmed(signal1$score[signal1$seqname == i], 9), type="h", bty="n", lwd=lwd, lend=2,
       xlab="", xaxt="n", ylim=ylim_LB1, col=cols[2], xlim=c(from,to), cex.axis=1.5, yaxt=yaxt_LB1, ylab=ylab_LB1, cex.lab=2, cex.axis=2)
  #mtext("Neuro-precursor cells", 3, -1)
  par(mar=c(0,5,15,5), new=TRUE)
  plot(0, type="n", bty="n", xlab="", ylab="", xaxt="n", yaxt="n", xlim=c(from, to), ylim=ylim_K36, yaxs="i")
  #rect(0, 0, 1e9, 100, col="blue", border=NA)
  lines(signal1_K36$start[signal1_K36$seqname == i], signal1_K36$score[signal1_K36$seqname == i], type="h", bty="n", lend=2, col="black")
  lines(c(ylim_K36[1], 1e9), rep(ylim_K36[1], 2), type="l", col=cols[2], lwd=3);
  if (axes_K36) {
    axis(4, cex.axis=2)
    mtext("H3K36me3", 4, line=3, cex=1.5)
  }
  
  par(mar=c(5,5,3,5))
  plot(signal2$start[signal2$seqname == i], runmed(signal2$score[signal2$seqname == i], 9), type="h", bty="n", lwd=lwd, lend=2,
       xlab="", xaxt="n", ylim=ylim_LB1, col=cols[3], xlim=c(from,to), cex.axis=1.5, yaxt=yaxt_LB1, ylab=ylab_LB1, cex.lab=2, cex.axis=2)
  #mtext("Astrocytes", 3, -1)
  
  #yaxt_LB1="s"; ylab_LB1=expression(log[2](LaminB1~interaction));

  #mtext(expression(log[2](LaminB1~interaction)), 2, -2, cex=1.5, at=0.55, outer=TRUE)
  mtext(ylab_text, 2, -2, cex=1.5, at=0.55, outer=TRUE)
  #mtext("H3K36me3 (ChIP-seq data)", 4, -2, cex=1.5, at=0.53, outer=TRUE)

  #par(mar=c(5,5,0,5))
  par(new=TRUE)
  par(mar=c(0,5,3,5))
  options("scipen"=5)
  plot(signal2$start[signal2$seqname == i], runmed(signal2$score[signal2$seqname == i], 49), 
       type="n", bty="n", ylim=c(0.4, 2.1), xlab="", ylab="", xlim=c(from,to), yaxt="n", cex.axis=2)
  options("scipen"=0)
  mtext(paste("chromosomal location (bp, ", chr, ")", sep=""), 1, 3.5, cex=1.5)

  rect(Entrez_full$start[Entrez_full$seqname == i & Entrez_full$strand == "+"], 0.77, Entrez_full$end[Entrez_full$seqname == i & Entrez_full$strand == "+"], 0.78, border=NA, col=1);
  rect(Entrez_full$start[Entrez_full$seqname == i & Entrez_full$strand == "-"], 0.57, Entrez_full$end[Entrez_full$seqname == i & Entrez_full$strand == "-"], 0.58, border=NA, col=1);

  text(Entrez_full$start[Entrez_full$seqname == i & Entrez_full$strand == "+"], 0.88, labels=Entrez_full$frame[Entrez_full$seqname == i & Entrez_full$strand == "+"], adj=c(0,0), col=1, cex=2)
  text(Entrez_full$end[Entrez_full$seqname == i & Entrez_full$strand == "-"], 0.68, labels=Entrez_full$frame[Entrez_full$seqname == i & Entrez_full$strand == "-"], adj=c(1,0), col=1, cex=2)
  
  rect(UTRs_full$start[UTRs_full$seqname == i & UTRs_full$strand == "+"], 0.74, UTRs_full$end[UTRs_full$seqname == i & UTRs_full$strand == "+"], 0.81, border=NA, col=1);
  rect(UTRs_full$start[UTRs_full$seqname == i & UTRs_full$strand == "-"], 0.54, UTRs_full$end[UTRs_full$seqname == i & UTRs_full$strand == "-"], 0.61, border=NA, col=1);
  
  rect(CDSs_full$start[CDSs_full$seqname == i & CDSs_full$strand == "+"], 0.7, CDSs_full$end[CDSs_full$seqname == i & CDSs_full$strand == "+"], 0.85, border=NA, col=1);
  rect(CDSs_full$start[CDSs_full$seqname == i & CDSs_full$strand == "-"], 0.5, CDSs_full$end[CDSs_full$seqname == i & CDSs_full$strand == "-"], 0.65, border=NA, col=1);

  # Again, only genes used in analysis, in red.
  if (emph_used_genes) {
    rect(Entrez$start[Entrez$seqname == i & Entrez$strand == "+"], 0.77, Entrez$end[Entrez$seqname == i & Entrez$strand == "+"], 0.78, border=NA, col=2);
    rect(Entrez$start[Entrez$seqname == i & Entrez$strand == "-"], 0.57, Entrez$end[Entrez$seqname == i & Entrez$strand == "-"], 0.58, border=NA, col=2);

    text(Entrez$start[Entrez$seqname == i & Entrez$strand == "+"], 0.88, labels=Entrez$frame[Entrez$seqname == i & Entrez$strand == "+"], adj=c(0,0), col=2, cex=2)
    text(Entrez$end[Entrez$seqname == i & Entrez$strand == "-"], 0.68, labels=Entrez$frame[Entrez$seqname == i & Entrez$strand == "-"], adj=c(1,0), col=2, cex=2)
  
    rect(UTRs$start[UTRs$seqname == i & UTRs$strand == "+"], 0.74, UTRs$end[UTRs$seqname == i & UTRs$strand == "+"], 0.81, border=NA, col=2);
    rect(UTRs$start[UTRs$seqname == i & UTRs$strand == "-"], 0.54, UTRs$end[UTRs$seqname == i & UTRs$strand == "-"], 0.61, border=NA, col=2);
  
    rect(CDSs$start[CDSs$seqname == i & CDSs$strand == "+"], 0.7, CDSs$end[CDSs$seqname == i & CDSs$strand == "+"], 0.85, border=NA, col=2);
    rect(CDSs$start[CDSs$seqname == i & CDSs$strand == "-"], 0.5, CDSs$end[CDSs$seqname == i & CDSs$strand == "-"], 0.65, border=NA, col=2);
  }

  #sel_genes <- Entrez_full$frame[Entrez_full$seqname == i & Entrez_full$start > from & Entrez_full$end < to]
  sel_genes <- Entrez$frame[Entrez$seqname == i & Entrez$start > from & Entrez$end < to]
  invisible(sel_genes)
}

# Used to shift doms half a inter-dom size up, i.e., rotate them circularly over the genome.
shiftDoms <- function(doms, seq_dims) {
  shift_doms <- NULL;
  for (chr in unique(doms$seqname)) {
    end <- seq_dims$end[seq_dims$seqname==chr]
    dom_chr <- doms[doms$seqname == chr,]
    n <- nrow(dom_chr)

    ctr <- floor((dom_chr$start + dom_chr$end)/2)
    ictr <- floor((dom_chr$end[1:(n-1)] + dom_chr$start[2:n])/2)
    ictr <- c(ictr, floor((dom_chr$end[n] + end)/2));

    dom_chr$start <- ctr;
    dom_chr$end <- ictr;

    shift_doms <- rbind(shift_doms, dom_chr);
  }
  invisible(shift_doms)
}

combineAlign <- function(aln1, aln2) {
  aln2$aln$relDomIdx <- aln2$aln$relDomIdx + nrow(aln1$relDom)
  aln <- list(aln    = rbind(aln1$aln, aln2$aln), 
              relDom = rbind(aln1$relDom, aln2$relDom));
  invisible(aln);
}

func.mem.use <- function() {
  sum_size <- 0;
  for (i in ls(envir=parent.frame())) {
    sum_size <- sum_size + object.size(eval(parse(text=i), envir=parent.frame()));
  }

  print(sum_size, units="auto")
  invisible(sum_size);
}

mem.use <- function() {
  sort(sapply(ls(envir=.GlobalEnv), function(x) object.size(eval(parse(text=x)))))
}

# Convert genomic locations back to probe indices
get_probe_coords <- function(GFF, template=NP_ES) {
  probes <- CGintvlOverlap(GFF, template, type="center", checkData=FALSE)
  ranges <- lapply(probes, function(x) lapply(x, range, na.rm=T))
  probe_coords <- NULL;
  for (i in names(ranges)) {
    probe_coords <- rbind(probe_coords, data.frame(start=unlist(lapply(ranges[[i]], "[[", 1)),
                                                   end=unlist(lapply(ranges[[i]], "[[", 2))))
  }

  probe_coords[probe_coords == Inf | probe_coords == -Inf] <- NA;

  invisible(probe_coords)
}

# Used to extend GFF based on number of array probes. Used in synteny analysis of intergenic regions.
ext_probe_coords <- function(probe_coords, template=NP_ES, upstr=25, downstr=25) {
  chrom_mins <- tapply(1:nrow(template), template$seqname, min)
  chrom_maxs <- tapply(1:nrow(template), template$seqname, max)

  probe_coords$start <- pmax(probe_coords$start - upstr,   
      chrom_mins[template$seqname[probe_coords$start]], na.rm=T)
  probe_coords$end   <- pmin(probe_coords$end   + downstr, 
      chrom_maxs[template$seqname[probe_coords$end]],   na.rm=T)

  invisible(probe_coords)
}

# Inverse of get_probe_coords()
get_chr_coords <- function(probe_coords, template=NP_ES) {
  res <- data.frame(seqname=template$seqname[probe_coords[,1]],
                    start=template$start[probe_coords[,1]],
                    end=template$end[probe_coords[,2]],
                    stringsAsFactors=FALSE)

  hits <- which(is.na(res$start))
  chrom_mins <- tapply(template$start, template$seqname, min)
#  res$start[hits] <- chrom_mins[res$seqname[hits]]

  hits <- which(is.na(res$end))
  chrom_maxs <- tapply(template$end, template$seqname, max)
#  res$end[hits] <- chrom_maxs[res$seqname[hits]]

  invisible(res)
}

get_dom_border_feats <- function(doms, feats, winsize=10000) {
  bndr <- CGintvlBndr(doms)
  bndrExt <- CGintvlExtend(bndr, upstr=winsize, downstr=0, chromStartEnd=chromSE_MM_NCBIm37)

  feats$center <- floor((feats$start + feats$end)/2)
  feats$start <- feats$end <- feats$center;

  # Return doms
  bndrfeat <- unlist(CGintvlOverlap(feats, bndrExt, type="minimal")) 

  # Create 0/1 vector, also taking care of duplicates
  res <- rep(FALSE, nrow(bndr))
  res[bndrfeat] <- TRUE;

  invisible(res)
}

get_dom_border_feats_outward <- function(doms, feats, winsize=10000) {
  bndr <- CGintvlBndr(doms)
  bndrExt <- CGintvlExtend(bndr, upstr=winsize, downstr=0, chromStartEnd=chromSE_MM_NCBIm37)

  feats5 <- feats[feats$strand=="-",]; feats5$start <- feats5$end # Select appropriate TSSs
  feats3 <- feats[feats$strand=="+",]; feats3$end <- feats3$start # Select appropriate TSSs
  
  bndr5Ext <- bndrExt[seq(1,nrow(bndrExt), 2),]
  bndr3Ext <- bndrExt[seq(2,nrow(bndrExt), 2),]

  # Return doms
  bndr5feat <- unlist(CGintvlOverlap(feats5, bndr5Ext, type="minimal"))
  bndr3feat <- unlist(CGintvlOverlap(feats3, bndr3Ext, type="minimal"))

  # Create 0/1 vector, also taking care of duplicates
  res <- rep(FALSE, nrow(bndrExt))
  res[seq(1, nrow(bndrExt), 2)][bndr5feat] <- TRUE
  res[seq(2, nrow(bndrExt), 2)][bndr3feat] <- TRUE

  invisible(res)
}

# OK
get_bindata <- function(aln, bins=seq(-25*1e3, 25*1e3, length.out=11)) {
  aln$aln$center <- (aln$aln$start + aln$aln$end)/2

  # Retain only features of which the center lies within a bin.
  aln$aln <- aln$aln[aln$aln$center < max(bins) & aln$aln$center > min(bins),] 

  # Do actual binning, making sure the bins extend to the extremes
  binning <- cut(c(range(bins), aln$aln$center), breaks=bins)
  binning <- binning[3:length(binning)] # Remove extremes
  aln$aln$centerbin <- binning # Assign features to bins

  # For each domain, obtain features and calculate median value per bin.
  print("Building binned dataset...")
  binned <- tapply(1:nrow(aln$aln), aln$aln$relDomIdx, function(x) {
    perc <- round((max(x) / nrow(aln$aln) * 100), 0)
    cat(paste("\r", perc, "%", sep=""))
    tapply(aln$aln$score[x], aln$aln$centerbin[x], median, na.rm=T)
  })
  cat("\n");

  # Convert result into bin (this could probably be done more elegantly)
  binmat <- matrix(NA, nrow=length(binned), ncol=length(bins)-1) 
  i <- 1;
  feh <- lapply(binned, function(x) {
    binmat[i,] <<- x; 
    i <<- i+1
  })
  rownames(binmat) <- names(binned)
  colnames(binmat) <- levels(binning);

  invisible(list(binmat=binmat, limits=range(bins)))
}

get_probedata <- function(doms_probes, signal, numprobes=10, direction="both") {
  binmat_probes_prime5 <- t(apply(doms_probes, 1, function(x) x[1] + c(-numprobes:numprobes)))
  binmat_probes_prime3 <- t(apply(doms_probes, 1, function(x) x[2] - c(-numprobes:numprobes)))

  binmat_prime5 <- t(apply(binmat_probes_prime5, 1, function(x) signal$score[x]))
  binmat_prime3 <- t(apply(binmat_probes_prime3, 1, function(x) signal$score[x]))

  if (direction == "both") {
    n <- nrow(doms_probes)*2
    binmat <- matrix(NA, ncol=numprobes*2+1, nrow=n)
    binmat[seq(1, n/2)*2-1,] <- binmat_prime5;
    binmat[seq(1, n/2)*2,]   <- binmat_prime3;
  } else {
    n <- nrow(doms_probes)
    if (direction == "prime5") {
      binmat <- binmat_prime5;
    } else {
      binmat <- binmat_prime3;
    }
  }

  rownames(binmat) <- 1:n
  colnames(binmat) <- c(-numprobes:numprobes)

  invisible(list(binmat=binmat, limits=c(-numprobes,numprobes)*1200))
}

doFisherGeneric <- function(set1, set2) {
  A <- sum(set1 & set2, na.rm=T)
  B <- sum(set1 & !set2, na.rm=T)
  C <- sum(!set1 & set2, na.rm=T)
  D <- sum(!set1 & !set2, na.rm=T)

  tab <- matrix(c(A,C,B,D), ncol=2, dimnames=list(c("in set1","not in set1"), c("in set2","not in set2")))
  print(tab)
  perc_set1 <- A/(A+B)*100;
  perc_set2 <- C/(C+D)*100;
  print(paste("Set 1:", perc_set1))
  print(paste("Set 2:", perc_set2))

  res_twosided <- fisher.test(tab, alt="two.sided")$p.value
  print(res_twosided)

  print(paste("Enrichment ratio:", perc_set1/perc_set2))
  res_enrich <- fisher.test(tab, alt="greater")$p.value
  print(res_enrich)

  print(paste("Depletion ratio:", perc_set2/perc_set1))
  res_deplet <- fisher.test(tab, alt="less")$p.value
  print(res_deplet)

  invisible(list(perc_set1=perc_set1, perc_set2=perc_set2, ratio=perc_set1/perc_set2, 
                 p=res_twosided, p_enrich=res_enrich, p_deplet=res_deplet))
}

doHyperGeneric <- function(set1, set2) {
  q <- length(which(set1 & set2))
  m <- length(which(set2))
  n <- length(which(!set2))
  k <- length(which(set1))

  print(paste("Enrichment ratio:", (q / k) / ((m-q) / (m+n-k))))
  res_enrich <- phyper(q-1, m, n, k, lower.tail=FALSE)
  print(res_enrich)

  print(paste("Depletion ratio:", ((m-q) / (m+n-k)) / (q / k)))
  res_deplet <- phyper(q, m, n, k, lower.tail=TRUE)
  print(res_deplet)

  invisible(list(enrich=res_enrich, deplet=res_deplet))
}

likelihood.gauss <- function (parameters, sample) {
  sample <- sample[!is.na(sample) & is.finite(sample)];
  if (length(parameters) != 5)
     stop ("parameters must be vector of length 5")
  prop = parameters[1]
  if (prop < 0) prop = 0
  if (prop > 1) prop = 1
  m1 = parameters[2]
  s1 = parameters[3]
  m2 = parameters[4]
  s2 = parameters[5]

  #print(parameters)
  return (-sum(log( prop*dnorm(sample, m1, s1) + (1-prop)*dnorm(sample, m2, s2) )))
}

# Not in use anymore
find_center_optim <- function(data, plot=FALSE) {
  dens <- density(data, na.rm=T)
  slopeSign <- diff(dens$y) > 0
  slopeSignChange <- diff(slopeSign) < 0
  h <- which(slopeSignChange == TRUE) + 1
  print(paste(length(h), "modi found."))

  ord <- order(dens$y[h], decreasing=TRUE)
  sel <- which(ord %in% c(1,2))
  center <- dens$x[mean(h[sel])]

  out <- optim(c(0.5, dens$x[h[sel[1]]], 1, dens$x[h[sel[2]]], 1), likelihood.gauss, sample = data)

  if (plot) {
    plot(dens)
    abline(v=dens$x[h[sel]], col=2)
    abline(v=center, col=3)
    curve(dnorm(x, mean=out$par[2], sd=out$par[3]), col=4, add=T)
    curve(dnorm(x, mean=out$par[4], sd=out$par[5]), col=4, add=T)
  }

  invisible(out);
}

# Not in use anymore
find_center_cent <- function(data, plot=FALSE) {
  dens <- density(data, na.rm=T)
  slopeSign <- diff(dens$y) > 0
  slopeSignChange <- diff(slopeSign) < 0
  h <- which(slopeSignChange == TRUE) + 1
  print(paste(length(h), "modi found."))

  ord <- order(dens$y[h], decreasing=TRUE)
  sel <- which(ord %in% c(1,2))
  center <- dens$x[mean(h[sel])]

  if (plot) {
    plot(dens)
    abline(v=dens$x[h[sel]], col=2)
    abline(v=center, col=3)
  }

  invisible(center);
}

# Not in use anymore
find_center_dip_max <- function(data, plot=FALSE) {
  dens <- density(data, na.rm=T)
  slopeSign <- diff(-dens$y) > 0
  slopeSignChange <- diff(slopeSign) < 0
  h <- which(slopeSignChange == TRUE) + 1
  print(paste(length(h), "dips found."))
  dip <- h[which.max(dens$y[h])]
  #dip <- h[1]; # Use the first dip found.

  if (plot) {
    plot(dens)
    abline(v=dens$x[dip], col=2)
  }

  invisible(dens$x[dip]);
}

# Currently used method
find_center_dip <- function(data, plot=FALSE) {
  dens <- density(data, na.rm=T)

  slopeSign <- diff(dens$y) > 0
  slopeSignChange <- diff(slopeSign) < 0
  modi <- which(slopeSignChange == TRUE) + 1
  print(paste(length(modi), "modi found."))

  ord <- order(dens$y[modi], decreasing=TRUE)
  sel <- sort(ord[1:2])
  
  slopeSign <- diff(-dens$y) > 0
  slopeSignChange <- diff(slopeSign) < 0
  dips <- which(slopeSignChange == TRUE) + 1
  print(paste(length(dips), "dips found."))

  # Select minimal dip inbetween two maximum modi
  dip_sel <- which(dips > modi[sel[1]] & dips < modi[sel[2]]);
  dip <- dips[dip_sel[which.min(dens$y[dips[dip_sel]])]]

  if (plot) {
    plot(dens)
    abline(v=dens$x[dip], col=2)
  }

  invisible(dens$x[dip]);
}

# Not in use anymore
find_center_clust <- function(data, plot=FALSE) {
  data <- data[!is.na(data) & is.finite(data)]
  clust <- kmeans(data, centers=2)
  rng <- sort(unlist(tapply(data, clust$cluster, range)));
  center <- mean(rng[2:3])
  dens <- density(data, na.rm=T)

  if (plot) {
    plot(dens)
    abline(v=clust$centers, col=2)
    abline(v=center, col=3)
  }

  invisible(center);
}

cond.norm <- function (x, conditions = NULL, invariantset = NULL, dat.log.scale = TRUE, method = c("quantile", "global")){
    method <- method[1]
    method.available = c("quantile", "global")
    method.use <- method.available[charmatch(method, method.available)]
    if (is.na(method.use)) {
        stop(paste("This method, '", as.character(method), "', is not (yet) implemented. Please use 'quantile' or 'global'.",
            sep = ""))
    }
    if (method.use == "quantile") {
        if (is.null(invariantset)) {
            x.dist <- apply(apply(x, 1, sort), 1, median)
            x.norm <- t(apply(x, 1, function(xx, y) y[rank(xx)],
                y = x.dist))
        } else {
            n.col <- ncol(x)
            x.rank <- t(apply(x, 1, rank))
            if (length(invariantset) == 1) {
                invset <- order(apply(x.rank, 2, mad, constant = 1))[1:invariantset]
            }
            else {
                invset <- invariantset
            }
            mins <- where.is.min(x.rank)
            maxs <- where.is.max(x.rank)
            if (is.logical(invset)) {
                invset[c(mins, maxs)] <- TRUE
            }
            else {
                invset <- unique(c(invset, mins, maxs))
                invset <- sort(invset)
            }
            inv.x <- x[, invset]
            inv.x.dist <- apply(apply(inv.x, 1, sort), 1, median)
            x.norm <- t(apply(x, 1, function(y, invset, inv.x.dist) approx(sort(rank(y)[invset]),
                inv.x.dist, rank(y))$y, invset = invset, inv.x.dist = inv.x.dist))

            if (!is.null(conditions)) {
                if (nrow(x) != length(conditions)) {
                  stop("The length of the vector 'conditions' should correspond to the number of rows of the data. NOTE that the genes are in the columns and the conditions in the rows.")
                }
                cond.lev <- unique(conditions)
                for (i in cond.lev) {
                  which.lev <- (conditions == i)
                  if (sum(which.lev) > 1) {
                    x.lev <- x.norm[which.lev, ]
                    x.lev.dist <- apply(apply(x.lev, 1, sort),
                      1, median)
                    x.lev.norm <- t(apply(x.lev, 1, function(xx,
                      y) y[rank(xx)], y = x.lev.dist))
                    x.norm[which.lev, ] <- x.lev.norm
                  }
                }
            }
        }
    }
    else {
        if (method.use == "global") {
            x.mdn <- median(as.vector(x))
            x.mad <- median(apply(x, 1, mad))
            x.norm <- t(apply(x, 1, function(y) (y - median(y))/mad(y)))
            x.norm <- x.mdn + x.norm * x.mad
        }
        else {
            warning("Check that the function is correct! The method specified could not be found.")
        }
    }
    x.norm
}

### MUCH simpler versions of the cond.norm() function above.
### Also, this separates the empirical distribution from the actual normalisation, which can be convenient.
get_quant_norm_dist <- function (x) {
  dist <- apply(apply(x, 1, sort), 1, median)
  invisible(dist)
}

do_quant_norm <- function (x, dist=NULL) {
  if (is.null(dist)) dist <- get_quant_norm_dist(x);

  x.norm <- t(apply(x, 1, function(xx, y) y[rank(xx)], y = dist))
  invisible(x.norm)
}


plotfile <- function(filename="Rplot", type="png", width=7, height=7, device="bitmap", warn=FALSE) {
    format <- switch(type,
                     pdf = "pdfwrite",
                     png = "png16m",
                     type);
    if (device == "bitmap") {
        bitmap(file=paste(filename, "_", id(), ".", type, sep=""),
               type=format, res=144, taa=4, width=width, height=height)
    } else {
        eval(parse(text=device))(file=paste(filename, "_", id(), ".", type, sep=""),
                                 res=144, width=width, height=height)
    }
    if (warn) print("Do not forget to use 'dev.off()' afterwards!")
}


get_delta_stats_orig_single <- function(delta, NS_mean, NS_sd, probes, p.adjust.method="BH", signif=0.05) {
  ## Obtain significantly dissociated genes or intergenic regions
  lengths <- (probes$end-probes$start+1)
  means <- apply(probes, 1, function(x) { 
    if (is.na(x[1])) NA else
      mean(delta[x[1]:x[2]], na.rm=T) 
  })
  p_as <- pnorm(means, mean=NS_mean, sd=NS_sd/sqrt(2*lengths), lower.tail=FALSE)
  p_dis <- pnorm(means, mean=NS_mean, sd=NS_sd/sqrt(2*lengths), lower.tail=TRUE)
  hits_as <- (p.adjust(p_as, method=p.adjust.method) < signif)
  hits_dis <- (p.adjust(p_dis, method=p.adjust.method) < signif)

  list(p_as=p_as, p_dis=p_dis, hits_as=hits_as, hits_dis=hits_dis)
}

get_delta_stats_orig <- function(delta, NS_mean, NS_sd, probes, probes_IG,
                            p.adjust.method="BH", signif=0.05) {
  genes_res <- get_delta_stats_orig_single(delta, NS_mean, NS_sd, probes, p.adjust.method, signif)
  IG_res <- get_delta_stats_orig_single(delta, NS_mean, NS_sd, probes_IG, p.adjust.method, signif)

  list(p_as=genes_res$p_as, p_dis=genes_res$p_dis, hits_as=genes_res$hits_as, hits_dis=genes_res$hits_dis, 
       IG_p_as=IG_res$p_as, IG_p_dis=IG_res$p_dis, IG_hits_as=IG_res$hits_as, IG_hits_dis=IG_res$hits_dis)
}

## Obtain significantly dissociated/associated genes/IGs
get_delta_stats_single <- function(delta, NS_mean, NS_sd, probes, 
                                   p.adjust.method="BH", signif=0.05, acfs=NULL) {
  lengths <- (probes$end-probes$start+1)
  means <- apply(probes, 1, function(x) { 
    if (is.na(x[1])) NA else mean(delta[x[1]:x[2]], na.rm=T) 
  })

  if (is.null(acfs))
    acfs <- acf(delta, lag.max=max(lengths, na.rm=T), plot=FALSE)

  NS_var <- NS_sd^2

  covs <- acfs$acf * NS_var # Scale correlation by SD based on null distribution
  #covs <- acfs$acf * var(delta) # Scale correlation by SD based on delta distribution

  kcovs <- c(0, sapply(2:max(lengths, na.rm=T), function(k) sum(covs[2:k] * (k-1):1)))
  sumvars <- ((NS_var)*lengths) + (2 * kcovs[lengths])
  #sumvars <- (var(delta)*lengths) + (2 * kcovs[lengths])

  z <- ((means*lengths) - (NS_mean*lengths)) / sqrt(sumvars)

  p_as <- pnorm(means*lengths, mean=NS_mean*lengths, sd=sqrt(sumvars), lower.tail=FALSE)
  hits_as <- (p.adjust(p_as, method=p.adjust.method) < signif)

  p_dis <- pnorm(means*lengths, mean=NS_mean*lengths, sd=sqrt(sumvars), lower.tail=TRUE)
  hits_dis <- (p.adjust(p_dis, method=p.adjust.method) < signif)

  # For a two-tailed test do both the lower-tailed and the upper-tailed test and 
  # double the P-value of the smaller of the two results.
  p <- pmin(p_as, p_dis)*2
  hits <- (p.adjust(p, method=p.adjust.method) < signif)

  list(z=z, p=p, p_as=p_as, p_dis=p_dis, hits=hits, hits_as=hits_as, hits_dis=hits_dis)
}


## Obtain significantly dissociated/associated genes/IGs
get_delta_stats_single_mean <- function(delta, NS_mean, NS_sd, probes, 
                                        p.adjust.method="BH", signif=0.05, acfs=NULL) {
  lengths <- (probes$end-probes$start+1)
  means <- apply(probes, 1, function(x) { 
    if (is.na(x[1])) NA else mean(delta[x[1]:x[2]], na.rm=T) 
  })

  if (is.null(acfs))
    acfs <- acf(delta, lag.max=max(lengths, na.rm=T), plot=FALSE)

  NS_var <- NS_sd^2

  covs <- acfs$acf * NS_var # Scale correlation by SD based on null distribution
  #covs <- acfs$acf * var(delta) # Scale correlation by SD based on delta distribution
}

#
#### 1
#  kcovs <- NS_var * 2 * c(0, sapply(2:max(lengths, na.rm=T), function(k) sum(covs[2:k] * (k-1):1)))
#  sumvars <- (lengths*NS_var) + kcovs[lengths]
#  p_as <- pnorm(means*lengths, mean=NS_mean*lengths, sd=sqrt(sumvars), lower.tail=FALSE)
#  p_as[1:10]
#
#
#### 2
#  kcovs <- NS_var * 2 * c(0, sapply(2:max(lengths, na.rm=T), function(k) sum(covs[2:k] * (k-1):1)))
#  sumvars <- (NS_var/lengths) + (kcovs[lengths]/(lengths^2))
#  p_as <- pnorm(means, mean=NS_mean, sd=sqrt(sumvars), lower.tail=FALSE)
#  p_as[1:10]
#
#
#
#### OK
#  kcovs <- c(0, sapply(2:max(lengths, na.rm=T), function(k) sum(covs[2:k] * (k-1):1)))
#  sumvars <- (NS_var*(lengths+2*kcovs[lengths]))
#  p_as <- pnorm(means*lengths, mean=NS_mean*lengths, sd=sqrt(sumvars), lower.tail=FALSE)
#  p_as[1:10]
#
#### OK
#  kcovs <- c(0, sapply(2:max(lengths, na.rm=T), function(k) sum(covs[2:k] * (k-1):1)))
#  sumvars <- (NS_var/lengths) + (NS_var*2*kcovs[lengths])/(lengths^2)
#  p_as <- pnorm(means, mean=NS_mean, sd=sqrt(sumvars), lower.tail=FALSE)
#  p_as[1:10]
#
#### OK
#  kcovs <- c(0, sapply(2:max(lengths, na.rm=T), function(k) sum(covs[2:k] * (k-1):1)))
#  sumvars <- (NS_var + (NS_var*2*kcovs[lengths])/lengths) / lengths
#  p_as <- pnorm(means, mean=NS_mean, sd=sqrt(sumvars), lower.tail=FALSE)
#  p_as[1:10]
#
#### OK
#  kcovs <- c(0, sapply(2:max(lengths, na.rm=T), function(k) sum(covs[2:k] * (k-1):1)))
#  sumvars <- (NS_var * (lengths + 2*kcovs[lengths]))/lengths^2
#  p_as <- pnorm(means, mean=NS_mean, sd=sqrt(sumvars), lower.tail=FALSE)
#  p_as[1:10]
#
#### OK
#  kcovs <- c(0, sapply(2:max(lengths, na.rm=T), function(k) sum(covs[2:k] * (k-1):1)/(k/2)))
#  sumvars <- (NS_var * (1 + kcovs[lengths]))/lengths
#  p_as <- pnorm(means, mean=NS_mean, sd=sqrt(sumvars), lower.tail=FALSE)
#  p_as[1:10]
#

get_delta_stats <- function(delta, NS_mean, NS_sd, probes, probes_IG,
                            p.adjust.method="BH", signif=0.05, acfs=NULL) {
  genes_res <- get_delta_stats_single(delta, NS_mean, NS_sd, probes, p.adjust.method, signif, acfs)
  IG_res <- get_delta_stats_single(delta, NS_mean, NS_sd, probes_IG, p.adjust.method, signif, acfs)

  ## 'Correct' NA results for intergenic regions: these result from the IGs being too small to calculate significance.
  #IG_res$hits[is.na(IG_res$hits)] <- TRUE
  #IG_res$hits_as[is.na(IG_res$hits_as)] <- TRUE
  #IG_res$hits_dis[is.na(IG_res$hits_dis)] <- TRUE

  list(z=genes_res$z, p=genes_res$p, p_as=genes_res$p_as, p_dis=genes_res$p_dis, 
       hits=genes_res$hits, hits_as=genes_res$hits_as, hits_dis=genes_res$hits_dis,
       IG_z=IG_res$z, IG_p=IG_res$p, IG_p_as=IG_res$p_as, IG_p_dis=IG_res$p_dis, 
       IG_hits=IG_res$hits, IG_hits_as=IG_res$hits_as, IG_hits_dis=IG_res$hits_dis)
}

## Obtain significantly dissociated/associated running windows
get_delta_stats_win <- function(delta, NS_mean, NS_sd, k=10, 
                                p.adjust.method="BH", signif=0.05) {
  means <- runmean(delta, k=k);
 
  acfs <- acf(delta, lag.max=k, plot=FALSE)

  NS_var <- NS_sd^2

  covs <- acfs$acf * NS_var # Scale correlation by SD based on null distribution
  #covs <- acfs$acf * var(delta) # Scale correlation by SD based on delta distribution

  kcovs <- c(0, sapply(2:k, function(l) sum(covs[2:l] * (l-1):1)))
  sumvars <- ((NS_var)*k) + (2 * kcovs[k])
  #sumvars <- (var(delta)*k) + (2 * kcovs[k])

  z <- ((means*k) - (NS_mean*k)) / sqrt(sumvars)

  p_as <- pnorm(means*k, mean=NS_mean*k, sd=sqrt(sumvars), lower.tail=FALSE)
  hits_as <- (p.adjust(p_as, method=p.adjust.method) < signif)

  p_dis <- pnorm(means*k, mean=NS_mean*k, sd=sqrt(sumvars), lower.tail=TRUE)
  hits_dis <- (p.adjust(p_dis, method=p.adjust.method) < signif)

  # For a two-tailed test do both the lower-tailed and the upper-tailed test and 
  # double the P-value of the smaller of the two results.
  p <- pmin(p_as, p_dis)*2
  hits <- (p.adjust(p, method=p.adjust.method) < signif)

  list(z=z, p=p, p_as=p_as, p_dis=p_dis, hits=hits, hits_as=hits_as, hits_dis=hits_dis)
}


Davies.Harte.FD <- function (n, sigma = 1, delta = 0) {
    N <- 2^ceiling(log2(n))
    M <- 2*N

    S.X.tau <- sigma^2 / pi * sin(pi*delta) * gamma(1-2*delta) *
        exp(lgamma((0:N)+delta) - lgamma((0:N)+1-delta))

    if (delta < 0)
        S.X.tau[1] <- -S.X.tau[1]

    S <- Re(fft(c(S.X.tau, S.X.tau[N:2])))

    Z <- rnorm(M)

    Y <- double(M)
        Y[1] <- Z[1] * sqrt(M*S[1])
        Y[2:N] <- (Z[2*(2:N)-1] + Z[2*(2:N)]*1i) * sqrt(M*S[2:N]/2)
        Y[N+1] <- Z[M] * sqrt(M*S[N+1])
        Y[(N+2):M] <- Conj(Y[N:2])

    Y <- Re(fft(Y, inverse = T)) / M

    return(Y[1:n])
} 

## Make synteny mapping between master (e.g., mouse) and slave (e.g., human) species.
make_synteny_maps <- function(master, slave, master_template, slave_template, 
                              master_species="mouse", slave_species="human", 
                              format="png", sample=TRUE, info="", master_col="black", 
                              slave_col=NULL, lwd=0.2, ylim=NULL, to_file=TRUE, outputdir="misc_figures", 
                              postfun=NULL, postpar=NULL, prefun=NULL, prepar=NULL) {

  master_chroms <- chrom_info[[master_species]];
  master_chromsizes <- chromsize_info[[master_species]];

  ncols <- round(sqrt(length(master_chroms)))
  nrows <- ceiling(sqrt(length(master_chroms)))
  tot <- ncols * nrows;
  colwidths <- tapply(master_chroms, rep(1:ncols, each=nrows)[1:length(master_chroms)], 
                      function(x) max(master_template$start[master_template$seqname %in% x]))

  if (to_file) {
    type <- format;
    if (format == "pdf") type <- "pdfwrite";
    if (format == "png") type <- "png16m";
    bitmap(file=paste(outputdir, "/synteny_maps_", info, "_", id(), ".", format, sep=""), 
           type=type, width=14*ncols, height=7*nrows, res=144, taa=4)
  }

  par(oma=c(2,3,2,0), mar=c(0,2,0,2))
  nf <- layout(matrix(1:(tot*3), 3*nrows, ncols, byrow = FALSE), heights=rep(c(3,3,1.5), tot), widths=colwidths/sum(colwidths))
  
  if (length(master_col) == 1) master_col <- rep(master_col, length(master))
  
  if (is.null(slave_col)) {
    slave_col <- rep("#000000", length(unique(slave_template$seqname)))
    names(slave_col) <- sort(unique(slave_template$seqname));
    slave_chroms <- chrom_info[[slave_species]];
    slave_col[slave_chroms] <- rainbow(length(slave_chroms));
    slave_col=slave_col[slave_template$seqname]
  }

  #xlim <- range(master_chromsizes[, c("start","end")])
  colwidths_chr <- rep(colwidths, each=nrows);
  names(colwidths_chr) <- master_chroms;

  for (chr in master_chroms) {
    print(chr)
    sel <- which(master_template$seqname == chr)
    if (sample) sel <- sample(sel, round(length(sel)/5))

    xlim <- c(0, colwidths_chr[chr]);
    if (!is.null(prefun)) prefun(prepar);
    plot(master_template$start[sel], master[sel], type="h", xlim=xlim, ylim=ylim,
         xaxt="n", bty="n", xlab="", ylab="", cex.axis=0.8, lwd=lwd, col=master_col[sel])
    if (!is.null(postfun)) postfun(postpar);

    mtext(paste("  ", chr), 3, 0, cex=1, line=-0.2, adj=0)
    if (!is.null(prefun)) prefun(prepar);
    plot(master_template$start[sel], slave[sel], type="h", xlim=xlim, ylim=ylim,
         xaxt="n", bty="n", xlab="", ylab="", cex.axis=0.8, lwd=lwd, col=slave_col[sel])
    if (!is.null(postfun)) postfun(postpar);
  
    axDat <- .prettyTicks(c(0, max(master_template$start[master_template$seqname == chr])), ntick = 7)
    plot(master_template$start[master_template$seqname == chr], 
         rep(1, length(which(master_template$seqname == chr))), type="n", 
         bty="n", axes=FALSE, ylim=c(-0.3, 1), xlab="", ylab="", xlim=xlim)
    axis(1, at = axDat$at, labels = axDat$labels, cex.axis=1, line=-2.5, cex.axis=0.8, padj=-1.2)
  }
  
  mtext(paste("genomic location (", axDat$xunit, ")", sep=""), 1, line=0.5, cex=1.5, outer=TRUE)
  mtext(text=expression(log[2](LaminB1~interaction)), 2, 1, cex=1.5, outer=TRUE)
  
  if (to_file) dev.off()
}

make_synteny_maps_chr <- function(master, slave, master_template, slave_template, 
                                  master_species="mouse", slave_species="human", 
                                  chr="chr1", xlim=NULL, ylim=NULL, format="png", 
                                  sample=TRUE, main="", info="", lwd=0.2,
                                  master_col="black", slave_col=NULL,
                                  width=14, to_file=TRUE, outputdir="misc_figures", 
                                  prefun=NULL, prepar=NULL, postfun=NULL, postpar=NULL) {
  master_chroms <- chrom_info[[master_species]];
  master_chromsizes <- chromsize_info[[master_species]];

  if (is.null(xlim)) {
    xrange <- "full";
    xlim <- as.numeric(master_chromsizes[master_chromsizes$seqname==chr,c("start","end")])
  } else {
    xrange <- paste(xlim, collapse="-")
  }

  if (is.null(ylim)) {
    ylim_master <- quantile(master, prob=c(0.001,0.999), na.rm=T)
    ylim_slave  <- quantile(slave, prob=c(0.001,0.999), na.rm=T)
  } else {
    ylim_master <- ylim_slave <- NULL;
  }

  if (to_file) {
    type <- format;
    if (format == "pdf") type <- "pdfwrite";
    if (format == "png") type <- "png16m";
    bitmap(file=paste(outputdir, "/synteny_maps_", chr, "_", xrange, "_", 
                      info, "_", id(), ".", format, sep=""), 
           type=type, width=width, height=3, res=144, taa=4)
  }

  nf <- layout(matrix(1:3, 3, 1, byrow = FALSE), heights=c(3,3,1.5))
  #par(oma=c(2,4,2,0), mar=c(0,2,0,0))
  par(mar=c(0,7,0,1), oma=c(0,0,ifelse(nchar(main) > 0, 2, 1),0))
  
  if (length(master_col) == 1) master_col <- rep(master_col, length(master))
  
  if (is.null(slave_col)) {
    slave_col <- rep("#000000", length(unique(slave_template$seqname)))
    names(slave_col) <- sort(unique(slave_template$seqname));
    slave_chroms <- chrom_info[[slave_species]];
    slave_col[slave_chroms] <- rainbow(length(slave_chroms));
    slave_col=slave_col[slave_template$seqname]
  }

  sel <- which(master_template$seqname == chr)
  if (sample) sel <- sample(sel, round(length(sel)/5))


  plot(0, type="n", xlim=xlim, ylim=ylim_master, xaxt="n", bty="n", xlab="", ylab="", cex.axis=2)
  if (!is.null(prefun)) prefun(prepar);
  lines(master_template$start[sel], master[sel], type="h", lwd=lwd, col=master_col[sel])
  #mtext(paste("     ", chr), 3, 0, cex=2, line=-0.2, adj=0)
  if (!is.null(postfun)) postfun(postpar);

  plot(0, type="n", xlim=xlim, ylim=ylim_slave, xaxt="n", bty="n", xlab="", ylab="", cex.axis=2)
  if (!is.null(prefun)) prefun(prepar);
  lines(master_template$start[sel], slave[sel], type="h", lwd=lwd, col=slave_col[sel])
  if (!is.null(postfun)) postfun(postpar);
  
  axDat <- .prettyTicks(xlim, ntick = round(width/2))
  plot(master_template$start[master_template$seqname == chr], 
       rep(1, length(which(master_template$seqname == chr))), type="n", 
       bty="n", axes=FALSE, ylim=c(0, 1), xlab="", ylab="", xlim=xlim)
  axis(1, at = axDat$at, labels = axDat$labels, cex.axis=2.5, line=-7, padj=0.5)

  mtext(text=expression(log[2](nuclear~lamina~interaction)), side=2, -3, cex=2, at=0.6, outer=TRUE)
  mtext(paste("location on ", master_species, " chromosome ", sub("chr", "", chr), 
              " (", axDat$xunit, ")", sep=""), side=3, -7, cex=2)

  mtext(main, side=3, cex=2, outer=TRUE)
  
  if (to_file) dev.off()
}

make_ext_synteny_maps_chr <- function(master, slave, master_template, slave_template, 
                                      master_species="mouse", slave_species="human", 
                                      chr="chr1", xlim=NULL, format="png", sample=TRUE, 
                                      info="", master_col="black", slave_col=NULL,
                                      lwd=0.2, ylim=NULL, width=14, to_file=TRUE,
                                      outputdir="misc_figures",
                                      ext_dat=NULL, postfun=NULL, postpar=NULL, 
                                      prefun=NULL, prepar=NULL) {
  master_chroms <- chrom_info[[master_species]];
  master_chromsizes <- chromsize_info[[master_species]];

  if (is.null(xlim)) {
    xrange <- "full";
    xlim <- as.numeric(master_chromsizes[master_chromsizes$seqname==chr,c("start","end")])
  } else {
    xrange <- paste(xlim, collapse="-")
  }

  if (is.null(ylim)) {
    ylim_master <- quantile(master, prob=c(0.001,0.999), na.rm=T)
    ylim_slave  <- quantile(slave, prob=c(0.001,0.999), na.rm=T)
    if (is.null(ext_dat)) {
      ylim_ext  <- NULL;
    } else {
      ylim_ext  <- lapply(ext_dat, function(x) quantile(x, prob=c(0.001,0.999), na.rm=T))
    }
  } else {
    ylim_master <- ylim_slave <- ylim_ext <- NULL;
  }

  if (to_file) {
    type <- format;
    if (format == "pdf") type <- "pdfwrite";
    if (format == "png") type <- "png16m";
    bitmap(file=paste(outputdir, "/synteny_maps_", chr, "_", xrange, "_", 
                      info, "_", id(), ".", format, sep=""), 
           type=type, width=width, height=3, res=144, taa=4)
  }

  nf <- layout(matrix(1:(length(ext_dat)+3), (length(ext_dat)+3), 1, byrow = FALSE), 
               heights=c(3,3,rep(1.5, length(ext_dat)),1.5))
  #par(oma=c(2,4,2,0), mar=c(0,2,0,0))
  par(mar=c(0,7,0,1), oma=c(0,0,1,0))
  
  if (length(master_col) == 1) master_col <- rep(master_col, length(master))
  
  if (is.null(slave_col)) {
    slave_col <- rep("#000000", length(unique(slave_template$seqname)))
    names(slave_col) <- sort(unique(slave_template$seqname));
    slave_chroms <- chrom_info[[slave_species]];
    slave_col[slave_chroms] <- rainbow(length(slave_chroms));
    slave_col=slave_col[slave_template$seqname]
  }

  sel <- which(master_template$seqname == chr)
  if (sample) sel <- sample(sel, round(length(sel)/5))

  plot(0, type="n", xlim=xlim, ylim=ylim_master, xaxt="n", bty="n", xlab="", ylab="", cex.axis=2)
  if (!is.null(prefun)) prefun(prepar);
  lines(master_template$start[sel], master[sel], type="h", lwd=lwd, col=master_col[sel])
  #mtext(paste("     ", chr), 3, 0, cex=2, line=-0.2, adj=0)
  legend("topleft", "(x,y)", inset=0.001, legend=master_species, cex=2, bty="n")
  if (!is.null(postfun)) postfun(postpar);

  plot(0, type="n", xlim=xlim, ylim=ylim_slave, xaxt="n", bty="n", xlab="", ylab="", cex.axis=2)
  if (!is.null(prefun)) prefun(prepar);
  lines(master_template$start[sel], slave[sel], type="h", lwd=lwd, col=slave_col[sel])
    legend("topleft", "(x,y)", inset=0.001, legend=slave_species, cex=2, bty="n")
  if (!is.null(postfun)) postfun(postpar);
  
  for (i in 1:length(ext_dat)) {
    plot(0, type="n", xlim=xlim, ylim=ylim_ext[[i]], xaxt="n", bty="n", xlab="", ylab="", cex.axis=2)
    if (!is.null(prefun)) prefun(prepar);
    lines(master_template$start[sel], ext_dat[[i]][sel], type="h", lwd=lwd, col="black")
    legend("topleft", "(x,y)", inset=0.001, legend=names(ext_dat)[i], cex=2, bty="n")
    if (!is.null(postfun)) postfun(postpar);
  }
  
  axDat <- .prettyTicks(xlim, ntick = round(width/2))
  plot(master_template$start[master_template$seqname == chr], 
       rep(1, length(which(master_template$seqname == chr))), type="n", 
       bty="n", axes=FALSE, ylim=c(0, 1), xlab="", ylab="", xlim=xlim)
  axis(1, at = axDat$at, labels = axDat$labels, cex.axis=2.5, line=-7, padj=0.5)

  mtext(paste("location on ", master_species, " chromosome ", sub("chr", "", chr), 
              " (", axDat$xunit, ")", sep=""), side=3, -7, cex=2)
  mtext(text=expression(log[2](nuclear~lamina~interaction)), 2, -3, cex=2, at=0.6, outer=TRUE)
  
  if (to_file) dev.off()
}


draw_rainbow_bar_outertext <- function(template, chroms, format="png", header="Human\nchroms", to_file=TRUE, 
                                       mar=c(2,2,2,6), info="", outputdir="misc_figures", main_cex=2.5) {
  if (to_file) {
    type <- format; 
    if (format == "pdf") type <- "pdfwrite";
    if (format == "png") type <- "png16m";
    bitmap(file=paste(outputdir, "/rainbow_bar_outertext_", info, "_", id(), ".", format, sep=""), 
           type=type, res=144, taa=4, height=7, width=1)
  }

  par(oma=c(0,0,4,0), mar=mar)

  chroms_cols <- rep("#000000", length(unique(template$seqname)))
  names(chroms_cols) <- sort(unique(template$seqname));
  chroms_cols[chroms] <- rainbow(length(chroms));

  n <- length(chroms);
  cols <- chroms_cols[chroms]

  coly1 <- rev((1:n) / n)
  coly2 <- coly1 - 1/n

  plot(0, type="n", xlim=c(0,1), ylim=c(0,1), xaxs="i", yaxs="i", axes=FALSE, xlab="", ylab="")

  rect(0, coly1, 1, coly2, col=cols, border=NA)
  grid(nx=1, ny=n, col="black", lty=1); box();
  mtext(side=3, header, cex=main_cex, font=2, adj=0.5, line=-0.75, outer=TRUE)

  tickMarks <- ((0:n)/n)+(1/(2*n)); tickMarks <- tickMarks[-length(tickMarks)]
  labels <- chroms
  axis(side=4, at=tickMarks, labels=rev(labels), las=2, cex.axis=2)
  #axis(side=1, at=tickMarks[sel], labels=c(2, 5, ">8"), cex.axis=2, padj=0.3)

  if (to_file) dev.off()
}

draw_rainbow_bar <- function(template, chroms, format="png", header="Human\nchroms", to_file=TRUE, 
                             mar=c(2,2,2,2), info="", outputdir="misc_figures", main_cex=2.5) {
  if (to_file) {
    type <- format; 
    if (format == "pdf") type <- "pdfwrite";
    if (format == "png") type <- "png16m";
    bitmap(file=paste(outputdir, "/rainbow_bar_", info, "_", id(), ".", format, sep=""), 
           type=type, res=144, taa=4, height=7, width=1)

    par(oma=c(0,0,4,0), mar=mar)
  }

  chroms_cols <- rep("#000000", length(unique(template$seqname)))
  names(chroms_cols) <- sort(unique(template$seqname));
  chroms_cols[chroms] <- rainbow(length(chroms));

  n <- length(chroms);
  cols <- chroms_cols[chroms]

  coly1 <- rev((1:n) / n)
  coly2 <- coly1 - 1/n

  plot(0, type="n", xlim=c(0,1), ylim=c(0,1), xaxs="i", yaxs="i", axes=FALSE, xlab="", ylab="")

  rect(0, coly1, 1, coly2, col=cols, border=NA)
  grid(nx=1, ny=n, col="black", lty=1); box();
  mtext(side=3, header, cex=main_cex, font=2, adj=0.5, line=-0.75, outer=TRUE)

  text(0.5, coly1-(1/(2*n)), col="white", sub("chr", "", chroms), cex=2.5, font=2, adj=0.5)

  if (to_file) dev.off()
}

draw_rainbow_box <- function(template, chroms, format="png", header="Human chromosomes", to_file=TRUE, 
                             mar=c(2,2,2,2), info="", outputdir="misc_figures", main_cex=2) {
  n <- length(chroms);
  ncols <- round(sqrt(n))
  nrows <- ceiling(sqrt(n))
  tot <- ncols * nrows;

  if (to_file) {
    type <- format; 
    if (format == "pdf") type <- "pdfwrite";
    if (format == "png") type <- "png16m";
    bitmap(file=paste(outputdir, "/rainbow_box_", info, "_", id(), ".", format, sep=""), 
           type=type, res=144, taa=4, height=nrows/2, width=ncols)
  }

  par(mar=mar)

  mat <- matrix(1:tot, nrow=nrows, ncol=ncols);
  mat <- mat[nrows:1,]

  chroms_cols <- rep("#FFFFFF", tot);
  names(chroms_cols) <- chroms;
  chroms_cols[chroms] <- rainbow(n);
  
  image(mat, col=rev(chroms_cols), axes=FALSE);
  grid(nx=nrows, ny=ncols, col="black", lty=1); box();
  mtext(side=3, header, cex=main_cex, font=2, adj=0, line=0.5)

  text(x=rep((0:(nrows-1))/(nrows-1), ncols), 
       y=rep(((ncols-1):0)/(ncols-1), each=nrows), labels=chroms, cex=3) 

  if (to_file) dev.off()
}


draw_HMMcall_box <- function(calls=c("LAD", "iLAD"), cols=c("red", "green"), textcols=c("black", "black"),
                             format="png", header="Human lamina status", to_file=TRUE, mar=c(2,2,2,2), 
                             main_cex=2, info="", outputdir="misc_figures") {

  if (to_file) {
    type <- format; 
    if (format == "pdf") type <- "pdfwrite";
    if (format == "png") type <- "png16m";
    bitmap(file=paste(outputdir, "/HMMcall_box_", info, "_", id(), ".", format, sep=""), 
           type=type, res=144, taa=4, height=1, width=length(calls))
  }

  par(mar=mar)

  n <- length(calls)
  mat <- matrix(1:n, nrow=n);

  image(mat, col=cols, axes=FALSE);
  grid(nx=n, ny=1, col="black", lty=1); box();
  mtext(side=3, header, cex=main_cex, font=2, adj=0, line=0.5)

  text((0:(n-1))/(n-1), rep(0, n), calls, cex=2, col=textcols) 

  if (to_file) dev.off()
}

#violin.plot <- function (x, x.scaling = 1, autoscale = 0.995, 
#    individual.scaling = FALSE, median.bar = "grey50", x.pos, cex.names=1, ...) {
## Author: Guillaume Filion.
## Date: June 7, 2011.
## x: a list/data.frame of values.
## x.scaling: an x scaling factor to apply to the densities.
## autoscale: the minimum proportion of the data to be displayed on the plot.
## ...: density, angle, border, col, and lty are passed to polygon().
#
#    if (autoscale < 0 || autoscale > 1)
#        stop ("'autoscale' outside the interval [0,1]")
#
#    m <- length(x)
#    if (missing(x.pos))
#        x.pos <- 1:m
#    densities <- lapply(X = x, FUN = density, na.rm = TRUE)
#
#    # Compute nice spacing parameters between plots.
#    y.sup <- sapply(X = densities, FUN = function (x) max(x$y))
#    x.min <- min(sapply(X = densities, FUN = function (x) min(x$x)))
#    x.max <- max(sapply(X = densities, FUN = function (x) max(x$x)))
#    if(individual.scaling) 
#        y.sup <- y.sup/.8
#    else
#        y.sup <- rep(max(y.sup), length(y.sup))
#
#    extrargs = list(...)
#    # The density, angle, border, col, and lty parameters are passed
#    # to polygon(). The rest is passed to plot().
#    pol.index <- match(names(extrargs),
#        c("density", "angle", "border", "col", "lty"))
#    plotargs <- list()
#    polargs <- list()
#    if (length(extrargs) > 0) {
#        for (i in 1:length(extrargs)) {
#            if (is.na(pol.index[i]))
#                plotargs[[names(extrargs)[i]]] <- extrargs[[i]]
#            else
#                polargs[[names(extrargs)[i]]] <- extrargs[[i]]
#        }
#    }
#    # Set the axis labels to null if not specified.
#    if (is.null(plotargs[["xlab"]]))
#        plotargs[["xlab"]] <- ""
#    if (is.null(plotargs[["ylab"]]))
#        plotargs[["ylab"]] <- ""
#
#    # Find a nice y zoom if ylim is not specified.
#    if (is.null(plotargs[["ylim"]])) {
#        ymin <- min(sapply(X = x, FUN = quantile, prob = (1 - autoscale)/2,
#            na.rm = TRUE))
#        ymax <- max(sapply(X = x, FUN = quantile, prob = 0.5 + autoscale/2,
#            na.rm = TRUE))
#        plotargs[["ylim"]] <- c(ymin, ymax)
#    }
#
#    # Overwrite type, xaxt, x, and y.
#    plotargs[["type"]] <- 'n'
#    plotargs[["xaxt"]] <- 'n'
#    plotargs[["x"]] <- c(min(x.pos)-.6*x.scaling, max(x.pos)+.6*x.scaling)
#    plotargs[["y"]] <- c(x.min, x.max)
#
#    # Plot the frame.
#    eval(as.call(c(plot, plotargs)))
#
#    # Manually recycle the parameters passed to polygon().
#    if (is.null(polargs[["col"]]))
#        polargs[["col"]] <- 'black'
#    for (i in 1:length(polargs))
#        polargs[[i]] <- rep(eval(polargs[[i]]), length.out = m)
#
#    # Set the default border parameter to the value of col.
#    if (is.null(polargs[["border"]]))
#        polargs[["border"]] <- polargs[["col"]]
#    # Set color to black if not specified.
#    if (is.null(polargs[["col"]]))
#        polargs[["col"]] <- "black"
#
#    # Plot the x-axis.
#    axis(side = 1, at = 1:m, labels = names(x), cex.axis=cex.names)
#
#    # Plot the violins.
#    polargs.i <- list()
#    for (i in 1:m) { 
#        polargs.i[["x"]] <- x.scaling*(c(densities[[i]]$y,
#            -rev(densities[[i]]$y))) / (2.1*y.sup[i]) + x.pos[i]
#        polargs.i[["y"]] <- c(densities[[i]]$x, rev(densities[[i]]$x))
#        for (a in names(polargs))
#            polargs.i[[a]] <- polargs[[a]][i]
#        eval(as.call(c(polygon, polargs.i)))
#
#        # Plot the median bar.
#        if (!is.na(median.bar)) {
#            medians <- sapply(X = x, FUN = median, na.rm = TRUE)
#            segments(x0 = x.pos[i], x1 = x.pos[i]-.4, y0 = medians[i],
#                y1 = medians[i], lwd = 2, lty = 1, col = median.bar)
#            segments(x0 = x.pos[i], x1 = x.pos[i]+.4, y0 = medians[i],
#                y1 = medians[i], lwd = 2, lty = 1, col = median.bar)
#        }
#    }
#}

violin.plot <- function (x, names, x.scaling = 1, autoscale = 0.995, 
    individual.scaling = FALSE, median.bar = "grey50", x.pos, ...) {
# Author: Guillaume Filion.
# Date: June 7, 2011.
# x: a list/data.frame of values.
# x.scaling: an x scaling factor to apply to the densities.
# autoscale: the minimum proportion of the data to be displayed on the plot.
# ...: density, angle, border, col, and lty are passed to polygon().

    if (autoscale < 0 || autoscale > 1)
        stop ("'autoscale' outside the interval [0,1]")

    m <- length(x)
    if (missing(x.pos))
        x.pos <- 1:m
    densities <- lapply(X = x, FUN = density, na.rm = TRUE)

    # Compute nice spacing parameters between plots.
    y.sup <- sapply(X = densities, FUN = function (x) max(x$y))
    x.min <- min(sapply(X = densities, FUN = function (x) min(x$x)))
    x.max <- max(sapply(X = densities, FUN = function (x) max(x$x)))
    if(individual.scaling) 
        y.sup <- y.sup/.8
    else
        y.sup <- rep(max(y.sup), length(y.sup))

    extrargs = list(...)
    # The density, angle, border, col, and lty parameters are passed
    # to polygon(). The rest is passed to plot().
    pol.index <- match(names(extrargs),
        c("density", "angle", "border", "col", "lty"))
    plotargs <- list()
    polargs <- list()
    if (length(extrargs) > 0) {
        for (i in 1:length(extrargs)) {
            if (is.na(pol.index[i]))
                plotargs[[names(extrargs)[i]]] <- extrargs[[i]]
            else
                polargs[[names(extrargs)[i]]] <- extrargs[[i]]
        }
    }
    # Set the axis labels to null if not specified.
    if (is.null(plotargs[["xlab"]]))
        plotargs[["xlab"]] <- ""
    if (is.null(plotargs[["ylab"]]))
        plotargs[["ylab"]] <- ""

    # Find a nice y zoom if ylim is not specified.
    if (is.null(plotargs[["ylim"]])) {
        ymin <- min(sapply(X = x, FUN = quantile, prob = (1 - autoscale)/2,
            na.rm = TRUE))
        ymax <- max(sapply(X = x, FUN = quantile, prob = 0.5 + autoscale/2,
            na.rm = TRUE))
        plotargs[["ylim"]] <- c(ymin, ymax)
    }

    # Overwrite type, xaxt, x, and y.
    plotargs[["type"]] <- 'n'
    plotargs[["xaxt"]] <- 'n'
    plotargs[["x"]] <- c(min(x.pos)-.6*x.scaling, max(x.pos)+.6*x.scaling)
    plotargs[["y"]] <- c(x.min, x.max)

    # Plot the frame.
    eval(as.call(c(plot, plotargs)))

    # Manually recycle the parameters passed to polygon().
    if (is.null(polargs[["col"]]))
        polargs[["col"]] <- 'black'
    for (i in 1:length(polargs))
        polargs[[i]] <- rep(eval(polargs[[i]]), length.out = m)

    # Set the default border parameter to the value of col.
    if (is.null(polargs[["border"]]))
        polargs[["border"]] <- polargs[["col"]]
    # Set color to black if not specified.
    if (is.null(polargs[["col"]]))
        polargs[["col"]] <- "black"

    # Plot the x-axis.
    axis(side = 1, at = 1:m, labels = names)

    # Plot the violins.
    polargs.i <- list()
    for (i in 1:m) { 
        polargs.i[["x"]] <- x.scaling*(c(densities[[i]]$y,
            -rev(densities[[i]]$y))) / (2.1*y.sup[i]) + x.pos[i]
        polargs.i[["y"]] <- c(densities[[i]]$x, rev(densities[[i]]$x))
        for (a in names(polargs))
            polargs.i[[a]] <- polargs[[a]][i]
        eval(as.call(c(polygon, polargs.i)))

        # Plot the median bar.
        if (!is.na(median.bar)) {
            medians <- sapply(X = x, FUN = median, na.rm = TRUE)
            segments(x0 = x.pos[i], x1 = x.pos[i]-.4, y0 = medians[i],
                y1 = medians[i], lwd = 2, lty = 1, col = median.bar)
            segments(x0 = x.pos[i], x1 = x.pos[i]+.4, y0 = medians[i],
                y1 = medians[i], lwd = 2, lty = 1, col = median.bar)
        }
    }
}

get_pref_stats_one_side <- function(hits, IG_hits) {
  n <- length(hits)*2
  all_hits <- rep(NA, n)
  all_hits[seq(1, n, 2)] <- IG_hits
  all_hits[seq(2, n, 2)] <- hits
  
  gene_ind <- rep(c(0,1), n/2)
 
  res <- doFisherGeneric(all_hits, gene_ind)

  invisible(res)
}

get_pref_stats <- function(stats) {
  print("Dissociation")
  pref_dis <- get_pref_stats_one_side(stats$hits_dis, stats$IG_hits_dis)
  print("Association")
  pref_as <- get_pref_stats_one_side(stats$hits_as, stats$IG_hits_as)

  invisible(list(pref_dis=pref_dis, pref_as=pref_as))
}

get_adj_stats_one_side <- function(hits, IG_hits, VERBOSE=TRUE) {
  # In case intergenic regions are too small to test for significance, 
  # we previously assigned 'NA' to their coordinates. Here, we change their
  # as/dis status to 'TRUE', in order to be able to glue adjacent dis/as genes together.
  IG_hits[is.na(IG_hits)] <- TRUE;

  n <- length(hits)*2
  all_hits <- rep(NA, n)
  all_hits[seq(1, n, 2)] <- IG_hits
  all_hits[seq(2, n, 2)] <- hits
  
  sel <- which(all_hits)
  adj <- data.frame(start=sel[c(1, which(diff(sel) != 1)+1)],
                    end=sel[c(which(diff(sel) != 1), length(sel))]);
 
  adj_full <- adj;

  # If the start or end position is odd, it is in IG. 
  # Because we're intereseted in genes only, we correct this by adding 1 to start and/or subtracting 1 from end.
  adj$start <- ifelse(adj$start %% 2 == 1, adj$start+1, adj$start)
  adj$end <- ifelse(adj$end %% 2 == 1, adj$end-1, adj$end)
  
  sel <- which(adj$start %% 2 == 0 & adj$end %% 2 == 0 & adj$start == adj$end)

  if (VERBOSE) {
    cat("% of genes that move in isolation: ");
    cat(length(sel) / length(which(hits)) * 100)
    cat("\n");
    cat("% of these genes that move without their intergenic regions: ");
    cat(length(which(adj_full$start[sel] %% 2 == 0 & adj_full$end[sel] %% 2 == 0)) / length(sel) * 100)
    cat("\n");
    cat("% of all genes that move in isolation and without their intergenic regions: ");
    cat(length(which(adj_full$start[sel] %% 2 == 0 & adj_full$end[sel] %% 2 == 0)) / length(which(hits)) * 100)
    cat("\n")
  }

  # Remove IG-only regions and translate indices back to gene indices
  adj <- adj[adj$end-adj$start >= 0,]
  adj$start <- adj$start/2
  adj$end <- adj$end/2
  
  ls <- (adj$end-adj$start+1)
  if (VERBOSE) print(table(ls)*as.numeric(names(table(ls)))/length(which(hits))*100)
  adj$length <- ls

  invisible(list(adj=adj, adj_full=adj_full))
}

get_adj_stats <- function(stats, VERBOSE=TRUE) {
  if (VERBOSE) print("Dissociating genes")
  adj_dis <- get_adj_stats_one_side(stats$hits_dis, stats$IG_hits_dis, VERBOSE)
  if (VERBOSE) print("Associating genes")
  adj_as <- get_adj_stats_one_side(stats$hits_as, stats$IG_hits_as, VERBOSE)

  invisible(list(adj_dis=adj_dis$adj, adj_as=adj_as$adj,
                 adj_dis_full=adj_dis$adj_full, adj_as_full=adj_as$adj_full))
}

get_flanks <- function(probes, n=10, inner=FALSE) {
  probes_5prime <- probes_3prime <- probes;

  if (inner) {
    probes_5prime$end   <- probes$start + n;
    probes_3prime$start <- probes$end   - n;
  } else {
    probes_5prime$end   <- probes$start;
    probes_5prime$start <- probes$start - n;
    probes_3prime$end   <- probes$end   + n;
    probes_3prime$start <- probes$end;
  }

  probes_flanks <- rbind(probes_5prime, probes_3prime)
  probes_flanks <- probes_flanks[order(probes_flanks$start),]

  # Make sure there are no overlapping flanks: OK
  m <- nrow(probes_flanks)
  res <- quantile(probes_flanks$start[2:m] - probes_flanks$end[1:(m-1)])
  if (any(res < 0)) {
    warning("Overlapping flanks");
    print(res);
  }

  invisible(probes_flanks)
}

### Used to calculate probe spacing in GFF
probe_spacing <- function(GFF) {
  summary(unlist(tapply(1:nrow(GFF), GFF$seqname, function(x) 
                   GFF$start[x[-1]] - GFF$end[x[-length(x)]])))
}


### Used to write diagnostics to the 'diag' object.
add_diag <- function(name, value) {
  if (length(value) == 1) {
    cat(paste(name, ": ", value, "\n", sep=""))
  } else {
    cat(paste(name, ":\n", sep=""))
    print(value)
  }
  diag[[name]] <<- value;
}

### Perform circular permutation using two binary sets: 
### Rotate the one and calculate overlap with the other.
### Then, compare to 'real' overlap.
get_circ_perm_p <- function(set1, set2, n=length(set1), do_sample=FALSE) {
  if (is.logical(set1)) set1 <- which(set1);
  if (is.logical(set2)) set2 <- which(set2);

  if (do_sample) {
    offsets <- sort(sample(seq.int(n-1), 1e4))
  } else {
    offsets <- seq.int(n-1)
  }
  rnd_overlap <- c();
  for (i in 1:length(offsets)) {
    offset <- offsets[i];
    if ((i/length(offsets)*100) %% 10 == 0)
      cat(paste(i/length(offsets)*100, "%\r", sep=""))

    rnd <- (((set2+offset)-1) %% n)+1
    rnd_overlap <- c(rnd_overlap,  length(intersect(set1,  rnd)))
  }
  cat("\n");

  ovl <- length(intersect(set1, set2))
  p_le <- length(which(rnd_overlap <= ovl)) / length(rnd_overlap)
  p_gr <- length(which(rnd_overlap >  ovl)) / length(rnd_overlap)
  p <- min(p_le, p_gr)*2

  res <- list(ovl=ovl, n=length(rnd_overlap), p_le=p_le, p_gr=p_gr, p=p, summary(rnd_overlap))
}

get_conc_circ_perm_p <- function(calls1, calls2, n=length(calls1), do_sample=FALSE) {
  if (do_sample) {
    offsets <- sort(sample(seq.int(n-1)+1, 1e4))
  } else {
    offsets <- seq.int(n-1)+1
  }

  rnd_overlap <- c();
  for (i in 1:length(offsets)) {
    offset <- offsets[i];
    if ((i/length(offsets)*100) %% 10 == 0)
      cat(paste(i/length(offsets)*100, "%\r", sep=""))

    rnd <- calls2[c(offset:n,1:(offset-1))]
    rnd_overlap <- c(rnd_overlap, sum(calls1 == rnd))
  }
  cat("\n");

  ovl <- sum(calls1 == calls2)

  p_le <- length(which(rnd_overlap <= ovl)) / length(rnd_overlap)
  p_gr <- length(which(rnd_overlap >  ovl)) / length(rnd_overlap)
  p <- min(p_le, p_gr)*2

  res <- list(ovl=ovl, n=length(rnd_overlap), p_le=p_le, p_gr=p_gr, p=p, sum=summary(rnd_overlap))
}

# define wideScreen function
wideScreen <- function(ncols) {
  if (missing(ncols)) 
    ncols <- as.integer(Sys.getenv("COLUMNS"))
  options(width=ncols);
}

#### Function to read in (large) .bed files and merge overlapping regions.
# NOTE: simply ignore the fact that UCSC coordinates are 0-based.
# We won't notice this on our probe-level resolution anyway.
readFixMergeBED <- function(file, purge=c(), verbose=FALSE) {
  alldat <- NULL;
  cnt <- 0;
  step <- 1e4;
  rcnt <- 0;
  cracks <- c();

  n <- as.numeric(system(paste("cat \"", file, "\" | wc -l", sep=''), intern=TRUE))
  print(paste("Total number of lines:", n));

  for (i in 1:ceiling(n/step)) {
    if (verbose)
      print(paste("Processing lines", (cnt*step)+1, "to", min(n, (cnt+1)*step)))
    else
      cat(".");

    toread <- min(step, n - (cnt*step));
    dat <- readBED(file, nlines=toread, skip=cnt*step, purge=purge);

    # Fix cases where start > end (this happens sometimes..)
    hits <- which(dat$start > dat$end)
    tmp <- dat$start[hits];
    dat$start[hits] <- dat$end[hits];
    dat$end[hits] <- tmp;

    # Sort data
    dat <- dat[order(dat$seqname, dat$start),]

    ## We extend genes by +/- 30bp, as this is half of our probe width.
    ## This is to facilitate further merging of annotations with small (< 60bp) gaps.
    #dat$start <- dat$start - 30;
    #dat$end   <- dat$end   + 30;

    # Merge overlapping entries.
    dat <- Merge(dat)
    left <- round((nrow(dat)/toread)*100, 2)
    if (verbose) 
      print(paste(left, "percent left after merging"))

    cracks <- c(cracks, rcnt + nrow(dat));
    alldat <- rbind(alldat, dat);
    rcnt <- nrow(alldat);

    cnt <- cnt+1;
  }
  if (!verbose) cat("\n")

  tmps <- list();
  idxs <- c();
  for (i in 1:(length(cracks)-1)) {
    crack <- cracks[i];                                   # Select crack to be processed
    idx <- crack+c(-100:100);                             # Expand crack
    idx <- idx[idx %in% 1:nrow(alldat)]                   # Correct for boundaries of alldat
    idxs <- c(idxs, idx);                                 # Remembers positions
    tmp <- alldat[idx,];                                  # Create temp object
    tmp <- tmp[order(tmp$seqname, tmp$start),]            # Sort object
    tmps <- rbind(tmps, Merge(tmp))                       # Merge and append to tmp collection
  }
  alldat <- alldat[-idxs,]                                # Remove crack (neighbour) positions from alldat
  alldat <- rbind(alldat, tmps)                           # Append new merged data
  alldat <- alldat[order(alldat$seqname, alldat$start),]  # Sort all data.

  print(paste(nrow(alldat), "(", round((nrow(alldat) / n) * 100, 2), "%) rows left after merging"), sep="")

  invisible(alldat);
}

### Function to calculate confidence intervals ###NEED TO CHECK!
get_conf_int <- function(x, dist="norm") {
  if (dist == "norm")
    qnorm(0.975) * (sd(x) / sqrt(length(x)))
  else
    qt(0.975, df=length(x)-1) * (sd(x) / sqrt(length(x)))
}

### Function to scale up/down a heatmap. Adapted from the original R heatmap() function.
heatmap.scale <-
function (x, Rowv = NULL, Colv = if (symm) "Rowv" else NULL, 
    distfun = dist, hclustfun = hclust, reorderfun = function(d, 
        w) reorder(d, w), do_reorder=TRUE, add.expr, symm = FALSE, revC = identical(Colv, 
        "Rowv"), scale = c("row", "column", "none"), na.rm = TRUE, 
    margins = c(5, 5), ColSideColors, RowSideColors, cexRow = 0.2 + 
        1/log10(nr), cexCol = 0.2 + 1/log10(nc), labRow = NULL, 
    labCol = NULL, main = NULL, xlab = NULL, ylab = NULL, keep.dendro = FALSE, 
    verbose = getOption("verbose"), width=4, height=4, cex.lab=NULL, ...) 
{
    scale <- if (symm && missing(scale)) 
        "none"
    else match.arg(scale)
    if (length(di <- dim(x)) != 2 || !is.numeric(x)) 
        stop("'x' must be a numeric matrix")
    nr <- di[1]
    nc <- di[2]
    if (nr <= 1 || nc <= 1) 
        stop("'x' must have at least 2 rows and 2 columns")
    if (!is.numeric(margins) || length(margins) != 2) 
        stop("'margins' must be a numeric vector of length 2")
    doRdend <- !identical(Rowv, NA)
    doCdend <- !identical(Colv, NA)
    if (is.null(Rowv)) 
        Rowv <- rowMeans(x, na.rm = na.rm)
    if (is.null(Colv)) 
        Colv <- colMeans(x, na.rm = na.rm)
    if (doRdend) {
        if (inherits(Rowv, "dendrogram")) 
            ddr <- Rowv
        else {
            hcr <- hclustfun(distfun(x))
            ddr <- as.dendrogram(hcr)
            if (!is.logical(Rowv) || Rowv) 
                ddr <- reorderfun(ddr, Rowv)
        }
        if (nr != length(rowInd <- order.dendrogram(ddr))) 
            stop("row dendrogram ordering gave index of wrong length")
    }
    else rowInd <- 1L:nr
    if (doCdend) {
        if (inherits(Colv, "dendrogram")) 
            ddc <- Colv
        else if (identical(Colv, "Rowv")) {
            if (nr != nc) 
                stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
            ddc <- ddr
        }
        else {
            hcc <- hclustfun(distfun(if (symm) 
                x
            else t(x)))
            ddc <- as.dendrogram(hcc)
            if (!is.logical(Colv) || Colv) 
                ddc <- reorderfun(ddc, Colv)
        }
        if (nc != length(colInd <- order.dendrogram(ddc))) 
            stop("column dendrogram ordering gave index of wrong length")
    }
    else colInd <- 1L:nc

    # Force no-reorder when neccessary.
    if (!do_reorder) {
      colInd <- 1L:nc;
      rowInd <- 1L:nr;
    }

    x <- x[rowInd, colInd]
    labRow <- if (is.null(labRow)) 
        if (is.null(rownames(x))) 
            (1L:nr)[rowInd]
        else rownames(x)
    else labRow[rowInd]
    labCol <- if (is.null(labCol)) 
        if (is.null(colnames(x))) 
            (1L:nc)[colInd]
        else colnames(x)
    else labCol[colInd]
    if (scale == "row") {
        x <- sweep(x, 1, rowMeans(x, na.rm = na.rm))
        sx <- apply(x, 1, sd, na.rm = na.rm)
        x <- sweep(x, 1, sx, "/")
    }
    else if (scale == "column") {
        x <- sweep(x, 2, colMeans(x, na.rm = na.rm))
        sx <- apply(x, 2, sd, na.rm = na.rm)
        x <- sweep(x, 2, sx, "/")
    }
    lmat <- rbind(c(NA, 3), 2:1)
    #lwid <- c(if (doRdend) 1 else 0.05, 4)
    lwid <- c(if (doRdend) 1 else 0.05, width)
    #lhei <- c((if (doCdend) 1 else 0.05) + if (!is.null(main)) 0.2 else 0,4)
    lhei <- c((if (doCdend) 1 else 0.05) + if (!is.null(main)) 0.2 else 0, height)
    if (!missing(ColSideColors)) {
        if (!is.character(ColSideColors) || length(ColSideColors) != 
            nc) 
            stop("'ColSideColors' must be a character vector of length ncol(x)")
        lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 1)
        lhei <- c(lhei[1], 0.2, lhei[2])
    }
    if (!missing(RowSideColors)) {
        if (!is.character(RowSideColors) || length(RowSideColors) != 
            nr) 
            stop("'RowSideColors' must be a character vector of length nrow(x)")
        lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 1), 
            1), lmat[, 2] + 1)
        lwid <- c(lwid[1], 0.2, lwid[2])
    }
    lmat[is.na(lmat)] <- 0
    if (verbose) {
        cat("layout: widths = ", lwid, ", heights = ", lhei, 
            "; lmat=\n")
        print(lmat)
    }
    op <- par(no.readonly = TRUE)
    on.exit(par(op))
    layout(lmat, widths = lwid, heights = lhei, respect = TRUE)
    if (!missing(RowSideColors)) {
        par(mar = c(margins[1], 0, 0, 0.5))
        image(rbind(1L:nr), col = RowSideColors[rowInd], axes = FALSE)
    }
    if (!missing(ColSideColors)) {
        par(mar = c(0.5, 0, 0, margins[2]))
        image(cbind(1L:nc), col = ColSideColors[colInd], axes = FALSE)
    }
    par(mar = c(margins[1], 0, 0, margins[2]))
    if (!symm || scale != "none") 
        x <- t(x)
    if (revC) {
        iy <- nr:1
        ddr <- rev(ddr)
        x <- x[, iy]
    }
    else iy <- 1L:nr
    image(1L:nc, 1L:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + 
        c(0, nr), axes = FALSE, xlab = "", ylab = "", ...)
    axis(1, 1L:nc, labels = labCol, las = 2, line = -0.5, tick = 0, 
        cex.axis = cexCol)
    if (!is.null(xlab)) 
        mtext(xlab, side = 1, line = margins[1] - 1.25, cex=cex.lab)
    axis(4, iy, labels = labRow, las = 2, line = -0.5, tick = 0, 
        cex.axis = cexRow)
    if (!is.null(ylab)) 
        mtext(ylab, side = 4, line = margins[2] - 1.25, cex=cex.lab)
    if (!missing(add.expr)) 
        eval(substitute(add.expr))
    par(mar = c(margins[1], 0, 0, 0))
    if (doRdend) 
        plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
    else frame()
    par(mar = c(0, 0, if (!is.null(main)) 1 else 0, margins[2]))
    if (doCdend) 
        plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
    else if (!is.null(main)) 
        frame()
    if (!is.null(main)) 
        title(main, cex.main = 1.5 * op[["cex.main"]])
    invisible(list(rowInd = rowInd, colInd = colInd, Rowv = if (keep.dendro && 
        doRdend) ddr, Colv = if (keep.dendro && doCdend) ddc))
}

### Obtain optimal clustering for a distance matrix.
get_optimal_ordering <- function(distmat) {
  library(cba)

  hc <- hclust(distmat, method="average");
  optim <- order.optimal(distmat, hc$merge);

  invisible(as.numeric(optim$order));
}

### Peak finder
# Taken from https://gist.github.com/jamiefolson/5831746
which.peaks <- function(x,partial=TRUE,decreasing=FALSE){
  if (decreasing){
    if (partial){
      which(diff(c(TRUE,diff(x)<=0,FALSE))>0)
    }else {
      which(diff(diff(x)<=0)>0)
    }    
  }else {
    if (partial){
      which(diff(c(TRUE,diff(x)>=0,FALSE))<0)
    }else {
      which(diff(diff(x)>=0)<0)
    }

  }
}


