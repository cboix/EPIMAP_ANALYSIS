#!/usr/bin/R

# --------------------------
# Text processing functions:
# --------------------------

#' Capitalize string
#' 
#' @param x string to capitalize
#' @export
capitalize <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1,1)), substring(s, 2),
      sep="", collapse=" ")
}

#' Split text by widths
#' 
#' @param x string to split
#' @param x text width
#' @export
split.text = function(x, width=90){
    x = as.character(x)
    sl = as.numeric(gregexpr(" ", x)[[1]])
    slat = c(1, sl[which(diff(sl %% width) < 0)], nchar(x))
    x2 = ""
    for (i in 2:length(slat)){
        x2 = paste0(x2, substr(x, slat[i-1], slat[i]), "\n")
    }
    x2 = gsub(" \n ", "\n", x2)
    x2 = gsub("\n$","",x2)
    return(x2)
}


#' Extract numbers from a string
#'
#' Only looks for integers with commas at the moment.
#' 
#' @param x string to search
#' @return all numbers found in the string
#' @export
munge.nos = function(x){
    locs = gregexpr("[0-9][0-9,]*",x)[[1]]
    ids = as.numeric(locs)
    lns = attributes(locs)$match.length
    # Extract:
    nos = sapply(1:length(ids), function(j){ 
                     num = substr(x, ids[j], ids[j]+ lns[j] - 1)
                     as.numeric(gsub(",","",num)) })
    return(nos)
}


#' Space labels out for 1d
#' 
#' @param xx positions of text
#' @param box.pad box size/padding for text
#' @param lim limits of space to use
#' @export
space.1d = function(xx, box.pad=0.02, lim=c(0,1)){
    # Identify centers of runs, shift things away from them:
    dx = round(diff(xx), 4)
    xl = dx < box.pad
    rl = rle(xl)
    rl$lengths
    rdf = data.frame(l=rl$l, v=rl$v, cs=cumsum(rl$l))
    rdf$cent = rdf$cs - (rdf$l - 1) / 2
    pivots = rdf$cent[rdf$v] + 1
    x2 = xx
    nx = length(xx)
    # Can't just do this, might perturb the order.
    for (pv in pivots){
        if (round(pv,0) == pv){
            bot = pv - 1
            top = pv + 1
        } else {
            mx = mean(dx[c(pv - .5, pv + .5)])
            dx[pv - .5] = mx - box.pad / 2
            if (pv + .5 < nx - 1){ dx[pv + .5] = mx + box.pad / 2}
            bot = pv - .5 - 1
            top = pv + .5 + 1
        }
        # Go through and fix ripples:
        for (i in rev(1:bot)){
            if (i >= 1){
                if (x2[i + 1] - x2[i] <= box.pad){
                    x2[i] = x2[i + 1] - box.pad
                } else { break }
            }
        } 
        for (i in top:(nx-1)){
            if (i > 1 & i <= nx){
                if (x2[i] - x2[i - 1] <= box.pad){
                    x2[i] = x2[i - 1] + box.pad
                } else { break }
            }
        }
    }
    x2[x2 > lim[2]] = lim[2]
    x2[x2 < lim[1]] = lim[1]
    # Fix after lim, from outside in:
    for (i in 1:round(nx/2)){
        if (x2[i + 1] - x2[i] <= box.pad){
            x2[i + 1] = x2[i] + box.pad
        } else { break }
    } 
    for (i in rev(round(nx/2):nx)){
        if (x2[i] - x2[i - 1] <= box.pad){
            x2[i - 1] = x2[i] - box.pad
        } else { break }
    } 
    return(x2)
}
