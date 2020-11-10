## drop.tip.R (2019-11-07)

##   Remove Tips in a Phylogenetic Tree

## Copyright 2003-2019 Emmanuel Paradis, 2017-2018 Klaus Schliep, 2018 Joseph Brown

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

## This file has been modified from its original by Daniel S. Caetano (2020-10-30)
## The modification has the objective of attending a more specific usage of the method.
## The original file for drop.tip is here: https://github.com/cran/ape/blob/master/R/drop.tip.R
## The original file for collapse.single is here: https://github.com/cran/ape/blob/master/R/collapse.singles.R
## Please cite the 'ape' package, and its authors, on any usage of this code.

collapse_singles_track_nodes <- function(tree, phy_edge){
    ## Function adapted from ape::collapse.single.
    ## edge_hist is a copy of the edge matrix with an additional column to track the history of the nodes.

    n <- length(tree$tip.label)
    tree <- reorder(tree) # this works now
    e1 <- tree$edge[, 1]
    e2 <- tree$edge[, 2]
    ## Name these vectors so we can track the operations later.
    names(e1) <- as.character( phy_edge[,3] )
    names(e2) <- as.character( phy_edge[,3] )

    ## The nodes that appear only once are single tip nodes.
    tab <- tabulate(e1)
    if (all(tab[-c(1:n)] > 1)) return(tree) # tips are zero

    if (is.null(tree$edge.length)) {
        root.edge <- FALSE
        wbl <- FALSE
    } else {
        wbl <- TRUE
        el <- tree$edge.length
    }

    ## start with the root node:
    ROOT <- n + 1L
    while (tab[ROOT] == 1) {
        ## This will only work if we have a single tip node at the root.
        i <- which(e1 == ROOT)
        ROOT <- e2[i]
        if (wbl) {
            el <- el[-i]
        }
        e1 <- e1[-i]
        e2 <- e2[-i]
    }

    singles <- which(tabulate(e1) == 1)

    if (length(singles) > 0) {
        ## This means we have work to do to fix the nodes.
        ii <- sort(match(singles, e1), decreasing = TRUE)
        jj <- match(e1[ii], e2)
        for (i in 1:length(singles)) {
            ## This is changing the values but not updating the names.
            ## But should it update the names? I don't know.
            e2[jj[i]] <- e2[ii[i]]
            names( e2[jj[i]] )
            names( e2[ii[i]] )
            if (wbl) el[jj[i]] <- el[jj[i]] + el[ii[i]]
        }
        e1 <- e1[-ii]
        e2 <- e2[-ii]
        if (wbl) el <- el[-ii]
    }
    Nnode <- length(e1) - n + 1L

    oldnodes <- unique(e1)
    if (!is.null(tree$node.label)){
        tree$node.label <- tree$node.label[oldnodes - n]
    }
    newNb <- integer(max(oldnodes))
    newNb[ROOT] <- n + 1L
    sndcol <- e2 > n
    e2[sndcol] <- newNb[e2[sndcol]] <- n + 2:Nnode
    e1 <- newNb[e1]
    tree$edge <- unname( cbind(e1, e2, deparse.level = 0) )

    ## Update the phy_edge object.
    phy_edge <- cbind(tree$edge, as.numeric(names(e2)), deparse.level = 0)

    tree$Nnode <- Nnode
    if (wbl) {
        tree$edge.length <- el
    }

    ## Return the updated tree and the phy_edge.
    return( list(phy = tree, phy_edge = phy_edge) )
}

#' Drop fossil tips and record matching nodes
#'
#' @param phy a phylogeny with fossil tips
#'
#' @return a list with the pruned phylogeny and a table with the matching nodes
#' @export
#' @importFrom ape is.rooted reorder.phylo is.ultrametric
#' @importFrom FossilSim prune.fossil.tips
drop_fossils_track_node <- function(phy){
    ## Drop fossil lineages and track the correspondence between the nodes.

    ## If ultrametric, then return the original phylogeny.
    if( is.ultrametric( phy ) ){
        warning( "phy is already ultrametric, returning the phylogeny." )
        return( list(pruned_phy = phy, node_table = NA) )
    }

    ## Add defaults.
    trim.internal <- TRUE
    subtree <- FALSE
    rooted <- is.rooted(phy)
    collapse.singles <- TRUE

    ## Keep copy of original tip.labels:
    og_tip_labels <- phy$tip.label

    ## Find the fossil tips to drop.
    ext_phy <- prune.fossil.tips(tree = phy)
    ## List of tips to drop:
    fossil_tips <- phy$tip.label[ !phy$tip.label %in% ext_phy$tip.label ]
    Ntip <- length(phy$tip.label)
    ## The tips to drop:
    tip <- which(phy$tip.label %in% fossil_tips)

    wbl <- !is.null(phy$edge.length)

    phy <- reorder.phylo(phy)
    NEWROOT <- ROOT <- Ntip + 1
    Nnode <- phy$Nnode
    Nedge <- dim(phy$edge)[1]
    edge1 <- phy$edge[, 1] # local copies
    edge2 <- phy$edge[, 2] #
    keep <- !logical(Nedge) # Create the vector.

    ## delete the terminal edges given by `tip':
    keep[match(tip, edge2)] <- FALSE

    if (trim.internal) {
        ints <- edge2 > Ntip
        ## delete the internal edges that do not have anymore
        ## descendants (ie, they are in the 2nd col of `edge' but
        ## not in the 1st one)
        repeat {
            sel <- !(edge2 %in% edge1[keep]) & ints & keep
            if (!sum(sel)) break
            keep[sel] <- FALSE
        }
    }

    ## drop the edges
    ## Keep the history of which edges do not exist anymore.
    ## Third column is a tag only.
    ## The flag is a copy of the original phy$edge[,1]. These are the internal nodes.
    phy_edge <- cbind(phy$edge, as.numeric(phy$edge[,1]))
    phy_edge <- phy_edge[keep, ] ## Now these two matrices need to mirror.
    phy$edge <- phy_edge[,1:2]

    if (wbl) phy$edge.length <- phy$edge.length[keep]

    ## find the new terminal edges (works whatever 'subtree' and 'trim.internal'):
    TERMS <- !(phy$edge[, 2] %in% phy$edge[, 1])

    ## get the old No. of the nodes and tips that become tips:
    oldNo.ofNewTips <- phy$edge[TERMS, 2]

    n <- length(oldNo.ofNewTips) # the new number of tips in the tree

    ## the tips may not be sorted in increasing order in the
    ## 2nd col of edge, so no need to reorder $tip.label
    ## Keep the same order of the tips, but renumber them to start from 1.
    phy_edge[TERMS, 2] <- rank(phy_edge[TERMS, 2]) ## Mirror the operation.
    phy$edge[TERMS, 2] <- rank(phy$edge[TERMS, 2])
    ## fix by Thomas Sibley (2017-10-28):
    ## Drop the tip labels:
    if (length(tip)) phy$tip.label <- phy$tip.label[-tip]

    phy$Nnode <- dim(phy$edge)[1] - n + 1L # update phy$Nnode

    ## The block below renumbers the nodes so that they conform
    ## to the "phylo" format
    ## Keep the history of the matching nodes!
    newNb <- integer(Ntip + Nnode)
    newNb[NEWROOT] <- n + 1L
    sndcol <- phy$edge[, 2] > n

    ## With the base tree
    newNb[sort(phy$edge[sndcol, 2])] <- (n + 2):(n + phy$Nnode)
    phy$edge[sndcol, 2] <- newNb[phy$edge[sndcol, 2]]
    phy$edge[, 1] <- newNb[phy$edge[, 1]]
    storage.mode(phy$edge) <- "integer"
    ## Mirror of the above block
    newNb[sort(phy_edge[sndcol, 2])] <- (n + 2):(n + phy$Nnode)
    phy_edge[sndcol, 2] <- newNb[phy_edge[sndcol, 2]]
    phy_edge[, 1] <- newNb[phy_edge[, 1]]
    storage.mode(phy_edge) <- "integer"

    if (!is.null(phy$node.label)){ # update node.label if needed
        phy$node.label <- phy$node.label[which(newNb > 0) - Ntip]
    }

    ## Here we need a function that can track the nodes.
    coll_list <- collapse_singles_track_nodes(tree = phy, phy_edge = phy_edge)

    ## Make the correspondence node matrix
    corr_mat <- unique( coll_list$phy_edge[,c(3,1)] )
    corr_dat <- data.frame(old_nodes = corr_mat[,1], new_nodes = corr_mat[,2])
    ## Make sure it is ordered in function of the new nodes.
    corr_dat <- corr_dat[order(corr_dat$new_nodes),]

    ## Create a matching table of tips to append.
    new_tip_labels <- coll_list$phy$tip.label
    old_tips <- sapply(new_tip_labels, function(x) which( x == og_tip_labels ))
    match_tips <- data.frame(old_nodes = unname(old_tips)
                             , new_nodes = 1:length(new_tip_labels))

    ## Make a complete table, including the corresponce among the tips.
    corr_dat <- rbind(corr_dat, match_tips)

    return( list(pruned_phy = coll_list$phy, node_table = corr_dat) )
}
