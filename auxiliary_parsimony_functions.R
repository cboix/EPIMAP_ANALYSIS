#!/usr/bin/R
library(parallel)

# Modified pratchet function:
# TODO: Improve speed of pars --> mclapply of optim.parsimony
mod.pratchet = function(data, start=NULL, method="fitch", maxit=1000, 
                        minit=10, k=10, trace=1, all=FALSE, rearrangements="SPR", 
                        perturbation="ratchet", multicore=TRUE, mc.cores=NULL,
                        savefile=NULL,...){
    print("[STATUS] Initializing")
    if (multicore && is.null(mc.cores)) { mc.cores <- detectCores() - 1 }
    search_history <- c(FALSE, FALSE)
    eps <- 1e-08
    trace <- trace - 1
    uniquetree <- function(trees){
        print("[STATUS] Calculating uniquetree")
        k <- 1
        res <- trees[[1]]
        result <- list()
        result[[1]] <- res
        k <- 2
        trees <- trees[-1]
        while (length(trees) > 0) {
            class(trees) <- "multiPhylo"
            rf <- suppressMessages(RF.dist(res, trees, FALSE))
            if (any(rf == 0)) { trees <- trees[-which(rf == 0)] }
            if (length(trees) > 0) {
                res <- trees[[1]]
                result[[k]] <- res
                k <- k + 1
                trees <- trees[-1]
            }
        }
        result
    }
    if (search_history[1]) { start_trees <- list() }
    if (search_history[2]) { search_trees <- list() }

    # Start with a tree based on hamming distance of the data
    tree <- NULL
    mp <- Inf
    if (perturbation != "random_addition") {
        if (is.null(start)) 
            start <- optim.parsimony(nj(dist.hamming(data)), 
                                     data, trace = trace, method = method, rearrangements = rearrangements, ...)
        tree <- start
        if (!is.binary(tree)) { tree <- multi2di(tree) }
        data <- subset(data, tree$tip.label)
        attr(tree, "pscore") <- parsimony(tree, data, method = method, 
                                          ...)
        mp <- attr(tree, "pscore")
        if (trace >= 0) { print(paste("Best pscore so far:", attr(tree, "pscore"))) }
    }
    FUN <- function(data, tree, method, rearrangements, ...) optim.parsimony(tree, 
                                                                             data = data, method = method, rearrangements = rearrangements, 
                                                                             ...)
    result <- list()
    result[[1]] <- tree
    # How to return trees:
    on.exit({
        if (!all) { result <- tree } else { class(result) <- "multiPhylo" }
        if (length(result) == 1) { result <- result[[1]] }
        return(result)
    })

    kmax <- 1
    nTips <- length(tree$tip.label)
    # For each iteration:
    for (i in seq_len(maxit)) {
        it.time = proc.time()
        print(paste("[STATUS] Iteration", i))
        if (perturbation == "ratchet") {
            print("[STATUS] Calculating bootstrap (parsimony ratchet)")
            ptm = proc.time()
            bstrees <- mod.bootstrap.phyDat(data, FUN, tree = tree, bs = 1, trace = trace,
                                            method = method, rearrangements = rearrangements, 
                                            multicore=multicore, mc.cores=mc.cores, ...)
            print(proc.time() - ptm)

            print("[STATUS] Applying optim.parsimony to each tree")
            ptm = proc.time()
            # trees <- lapply(bstrees, optim.parsimony, data, trace = trace, 
            #                 method = method, rearrangements = rearrangements, ...)
            trees <- mclapply(bstrees, optim.parsimony, data, trace = trace, 
                            method = method, rearrangements = rearrangements, ...., mc.cores=mc.cores)
            # trees = foreach(i=1:length(bstrees)) %dopar% {
            #     optim.parsimony(bstrees[[i]], data, trace = trace, 
            #                     method = method, rearrangements = rearrangements, ....)
            # }
            print(proc.time() - ptm)

            if (search_history[1]) { start_trees[[i]] <- bstrees[[1]] }
            if (search_history[2]) { search_trees[[i]] <- trees[[1]] }
        }

        if (perturbation == "stochastic") {
            print("[STATUS] Performing stochastic switches")
            treeNNI <- rNNI(tree, floor(nTips/2))
            print("[STATUS] Applying optim.parsimony to each tree")
            trees <- optim.parsimony(treeNNI, data, trace = trace, 
                                     method = method, rearrangements = rearrangements, ...)
            trees <- list(trees)
            if (search_history[1]) { start_trees[[i]] <- treeNNI }
            if (search_history[2]) { search_trees[[i]] <- trees[[1]] }
        }

        if (perturbation == "random_addition") {
            print("[STATUS] Performing random addition")
            treeRA <- random.addition(data)
            print("[STATUS] Applying optim.parsimony to each tree")
            trees <- optim.parsimony(treeRA, data, trace = trace, 
                                     method = method, rearrangements = rearrangements, ...)
            trees <- list(trees)
            if (search_history[1]) { start_trees[[i]] <- treeRA }
            if (search_history[2]) { search_trees[[i]] <- trees[[1]] }
        }

        if (inherits(result, "phylo")) { m <- 1 } else { m <- length(result) }
        if (m > 0) { trees[2:(1 + m)] <- result[1:m] }

        print("[STATUS] Looking at max pscores")
        pscores <- sapply(trees, function(data) attr(data, "pscore"))
        mp1 <- min(pscores)
        # Count streak:
        if ((mp1 + eps) < mp) { kmax <- 1 } else { kmax <- kmax + 1 }
        mp <- mp1
        if (trace >= 0) { print(paste("Best pscore so far:", mp)) }
        ind <- which(pscores < mp + eps)
        if (length(ind) == 1) {
            result <- trees[ind]
            tree <- result[[1]]
            if (!is.null(savefile)){
                print("[STATUS] Saving tree")
                write.tree(tree, file=savefile)
            }
        } else {
            result <- uniquetree(trees[ind])
            l <- length(result)
            tree <- result[[sample(l, 1)]]
            if (!is.null(savefile)){
                print("[STATUS] Saving tree")
                write.tree(tree, file=savefile)
            }
        }
        print("[STATUS] Iteration time:")
        print(proc.time() - it.time)
        # Quit if streaks is at least k long and past min # iterations:
        if ((kmax >= k) && (i >= minit)) { (break)() }
    }
}




mod.bootstrap.phyDat = function (x, FUN, bs = 100, multicore = FALSE, 
                                 mc.cores = NULL, jumble = TRUE, ...){
    if (multicore && is.null(mc.cores)) { mc.cores <- detectCores() - 1 }
    weight <- attr(x, "weight")
    v <- rep(seq_along(weight), weight)
    BS <- vector("list", bs)
    for (i in 1:bs) BS[[i]] <- tabulate(sample(v, replace = TRUE), length(weight))

    if (jumble) {
        J <- vector("list", bs)
        l <- length(x)
        for (i in 1:bs){ J[[i]] <- list(BS[[i]], sample(l)) }
    }

    fitPar <- function(weights, data, ...) {
        ind <- which(weights > 0)
        data <- phangorn:::getRows(data, ind)
        attr(data, "weight") <- weights[ind]
        FUN(data, ...) }

    fitParJumble <- function(J, data, ...) {
        ind <- which(J[[1]] > 0)
        data <- phangorn:::getRows(data, ind)
        attr(data, "weight") <- J[[1]][ind]
        data <- subset(data, J[[2]])
        FUN(data, ...) }

    if (multicore) {
        print("Multicore bootstrap")
        if (jumble) { 
            res <- mclapply(J, fitParJumble, x, ..., mc.cores = mc.cores)
        } else {
            res <- mclapply(BS, fitPar, x, ..., mc.cores = mc.cores)
        }
    } else {
        print("Single-core bootstrap")
        if (jumble) {
            res <- lapply(J, fitParJumble, x, ...)
        } else {
            res <- lapply(BS, fitPar, x, ...)
        }
    }

    if (class(res[[1]]) == "phylo") {
        class(res) <- "multiPhylo"
        res <- .compressTipLabel(res)
    }
    res
}




mod.optim.fitch = function(tree, data, trace = 1, rearrangements = "SPR", ...){
    # Check data:
    if (!inherits(tree, "phylo")) { stop("tree must be of class phylo") }
    if (!is.binary(tree)) {
        tree <- multi2di(tree)
        attr(tree, "order") <- NULL
    }
    if (is.rooted(tree)) {
        tree <- unroot(tree)
        attr(tree, "order") <- NULL
    }
    if (is.null(attr(tree, "order")) || attr(tree, "order") == "cladewise") {
        tree <- reorder(tree, "postorder") }
    if (class(data)[1] != "phyDat") { stop("data must be of class phyDat") }

    rt <- FALSE
    dup_list <- NULL
    addTaxa <- FALSE
    tmp <- TRUE
    star_tree <- FALSE
    while (tmp) {
        nam <- names(data)
        data <- phangorn:::removeParsUninfoSites(data)
        p0 <- attr(data, "p0")
        if (attr(data, "nr") == 0) {
            star_tree <- TRUE
            (break)()
            tmp <- FALSE
        }
        dup <- map_duplicates(data)
        if (!is.null(dup)) {
            tree <- drop.tip(tree, dup[, 1])
            if (length(tree$tip.label) > 2) 
                tree <- unroot(tree)
            tree <- reorder(tree, "postorder")
            dup_list <- c(list(dup), dup_list)
            addTaxa <- TRUE
            data <- subset(data, setdiff(names(data), dup[, 1]))
        } else { (break)() }
    }

    nr <- attr(data, "nr")
    nTips <- as.integer(length(tree$tip.label))
    if (nTips < 5) 
        rearrangements <- "NNI"
    data <- subset(data, tree$tip.label, order(attr(data, "weight"), 
        decreasing = TRUE))
    dat <- phangorn:::prepareDataFitch(data)
    weight <- attr(data, "weight")
    m <- nr * (2L * nTips - 2L)
    on.exit({
        .C("fitch_free")
        if (addTaxa) {
            if (rt) tree <- ptree(tree, data)
            for (i in seq_along(dup_list)) {
                dup <- dup_list[[i]]
                tree <- add.tips(tree, dup[, 1], dup[, 2])
            }
            tree
        }
        if (length(tree$tip.label) > 2) tree <- unroot(tree)
        attr(tree, "pscore") <- pscore + p0
        return(tree)
    })
    .C("fitch_init", as.integer(dat), as.integer(nTips * nr), 
        as.integer(m), as.double(weight), as.integer(nr))
    tree$edge.length <- NULL
    swap <- 0
    iter <- TRUE
    if (nTips < 4) { iter <- FALSE }
    pscore <- phangorn:::fast.fitch(tree, nr)

    while (iter) {
        res <- phangorn:::fitch.nni(tree, dat)
        tree <- res$tree
        if (trace > 1) {
            cat("optimize topology: ", pscore + p0, 
                "-->", res$pscore + p0, "\n") }
        pscore <- res$pscore
        swap <- swap + res$swap
        if (res$swap == 0) {
            if (rearrangements == "SPR") {
                tree <- phangorn:::fitch.spr(tree, dat)
                psc <- phangorn:::fast.fitch(tree, nr)
                if (trace > 1) 
                  cat("optimize topology (SPR): ", pscore + p0, 
                    "-->", psc + p0, "\n")
                if (pscore < psc + 1e-06) 
                  iter <- FALSE
                pscore <- psc
            }
            else iter <- FALSE
        }
    }

    if (trace > 0) { 
        cat("Final p-score", pscore + p0, "after ", 
            swap, "nni operations \n") }
}

