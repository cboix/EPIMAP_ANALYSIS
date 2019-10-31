##########################################################################################################
#                                                                                                        #
# Method for generating a master list / Index of DNaseI hypersensitivity sites.                          #
# All code, implementation details and design by Wouter Meuleman and Eric Rynes.                         #
#                                                                                                        #
# Version: WM20180313                                                                                    #
#                                                                                                        #
##########################################################################################################


merge_overlap <- function(DHSs, peaks, type="summit", force=FALSE) {
  num <- 0;
  #core_starts <- pmax(DHSs$start, DHSs$summit+DHSs$conf95_low, na.rm=T)
  #core_ends <- pmin(DHSs$end, DHSs$summit+DHSs$conf95_high, na.rm=T)
  core_starts <- pmax(DHSs$start, DHSs$disp_low, na.rm=T)
  core_ends <- pmin(DHSs$end, DHSs$disp_high, na.rm=T)
  score_ord <- order(DHSs$score);
  for (i in score_ord) {
    idx <- c()
    if (type=="summit") {
      idx <- which(DHSs$start <= DHSs$summit[i] & DHSs$end >= DHSs$summit[i])
    } else if (type=="core") {
      idx <- which(core_starts <= core_ends[i] & core_ends >= core_starts[i])
    } else if (type=="any") {
      idx <- which(DHSs$start <= DHSs$end[i] & DHSs$end >= DHSs$start[i])
    }
    idx <- idx[which.max(DHSs$score[idx])]

    if (length(idx) == 0) { 
      print(paste("Summit misalignment for DHS on line", i, "-- skipping merge")); 
      next 
    }

    # If merging these DHSs results in lots of multi-peak contributions from samples... don't merge.
    tab <- tabulate(table(peaks$sampleID[which(peaks$ID %in% DHSs$ID[c(i,idx)])]))
    if (which.max(tab * 1:length(tab)) != 1 & sum(tab) > 10 & !force) {
      print(paste("Not merging DHSs", DHSs$ID[i], "and", DHSs$ID[idx], 
              "because of too many multi-peak sample contributions!"))
      next; 
    }

    if (idx != i) {
      peaks$ID[which(peaks$ID==DHSs$ID[i])] <- DHSs$ID[idx]
      DHSs[idx,] <- get_FWHM(peaks[which(peaks$ID==DHSs$ID[idx]),])
      DHSs[i,] <- rep(NA, ncol(DHSs))
      num <- num+1;
    }
  } 
  DHSs <- DHSs[which(!is.na(DHSs$seqname)),]
  DHSs <- DHSs[order(DHSs$seqname, DHSs$start, DHSs$end),]

  list(DHSs=DHSs, peaks=peaks, num=num);
}

get_min_max <- function(x, fun="max", pos="center") { 
  idx <- which(x == eval(parse(text=fun))(x, na.rm=T)); 
  if (pos=="center") {
    return(idx[ceiling(length(idx)/2)])
  } else if (pos=="left") {
    return(idx[1])
  } else if (pos=="right") {
    return(idx[length(idx)])
  }
}

get_FWHM <- function(res_peaks) {
  # Obtain range of values and tabulate to create pile-up
  xlim <- range(c(res_peaks$start, res_peaks$end))
  tab <- tabulate(unlist(apply(res_peaks, 1, function(x) x[2]:x[3]))-xlim[1]+1)
  tab_max <- get_min_max(tab, fun="max") # Identify maximal point in pile-up

  # Determine upstream and downstream boundaries
  us <- c(1:tab_max)[get_min_max(abs(tab[1:tab_max] - (max(tab)/2)), fun="min", pos="left")]
  ds <- c(tab_max:length(tab))[get_min_max(abs(tab[tab_max:length(tab)] - (max(tab)/2)), fun="min", pos="right")]

  start <- us+xlim[1]-1
  end <- ds+xlim[1]-1
  score <- sum(tapply(res_peaks$score, res_peaks$sampleID, max))
  numsamples <- length(unique(res_peaks$sampleID))
  numpeaks <- length(res_peaks$sampleID)
  summit <- round(median(res_peaks$wavelet_summit));
  #dispersion <- mad(res_peaks$wavelet_summit) # constant = 1.4826
  #dispersion <- mad(res_peaks$wavelet_summit, constant=1) # constant = 1
  #conf95 <- sort(res_peaks$wavelet_summit-summit)[pmax(1, qbinom(c(.025,.975), length(res_peaks$wavelet_summit), 0.5))]
  #conf95_low <- conf95[1];
  #conf95_high <- conf95[2];
  #disp <- range(res_peaks$wavelet_summit);
  disp <- quantile(res_peaks$wavelet_summit, prob=c(0.025, 0.975));
  disp_low <- max(min(res_peaks$wavelet_summit), round(disp[1]));
  disp_high <- min(max(res_peaks$wavelet_summit), round(disp[2]));

  # If the (wavelet) summit is outside the delineated region, we use the median center instead.
  if (summit < start | summit > end) {
    summit <- round(median((start+end)/2))
    #dispersion <- NA
    #conf95_low <- conf95_high <- NA
    disp_low <- disp_high <- NA
  }

  # Return resulting region with summary statistics
  data.frame(seqname=res_peaks$seqname[1], start=start, end=end, ID=res_peaks$ID[1], score=score, numsamples=numsamples, 
             numpeaks=numpeaks, width=end-start, summit=summit, disp_low=disp_low, disp_high=disp_high, stringsAsFactors=FALSE)
}

get_cut_points <- function(summit_loc) {
  xlim <- range(summit_loc);
  n <- diff(xlim)+1
  if (n < 20) return(c(0, n))

  # turn ~20bp resolution data into 20bp-window smoothed versions
  cont_dat <- tabulate(summit_loc)
  smooth_dat <- runmean(cont_dat, 20, endrule="constant")

  # Determine potential points to chop up peak clump
  #cutlevel <- max(smooth_dat)/2 # Threshold for minimal signal to consider for cutting
  #cutlevel <- min(1, median(smooth_dat)) # WM20180303
  cutlevel <- median(smooth_dat) # WM20180302b -- works well!
  idx <- which(diff(smooth_dat < cutlevel) != 0)
  idx <- c(0, idx, n)

  cuts <- c(0)
  if (length(idx) > 2 & max(smooth_dat) >= 1) { # Only determine cut-sites if the maximum signal observed is at least 1
    cuts <- c(cuts, sapply(1:(length(idx)-1), function(i) {
      min_pos <- get_min_max(smooth_dat[idx[i]:idx[i+1]], fun="min", pos="center")
      ifelse(min_pos == 1, NA, idx[i] + min_pos -1)
    }))
  }
  cuts <- cuts[!is.na(cuts)]
  cuts <- c(cuts, n);

  # Clean up near the ends
  if (sum(smooth_dat[1:cuts[2]]) < sum(smooth_dat)*0.01)
    cuts <- cuts[-2];
  if (sum(smooth_dat[cuts[length(cuts)-1]:cuts[length(cuts)]]) < sum(smooth_dat)*0.01)
    cuts <- cuts[-(length(cuts)-1)];

  if (is.unsorted(cuts)) stop("Unordered cut-point list -- this should not happen!")

  return(cuts)
}

