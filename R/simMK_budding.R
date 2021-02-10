#' Simulate discrete trait evolution with budding and exponential decay function.
#'
#' This function adds a probability of budding speciation and allows only for a negative relationship between rates of trait evolution and lineage-age. Function also allows for cladogenetic change in the trait (a deterministic change in the state of the trait after budding speciation happens).
#' \cr
#' MOTHER LINEAGES: The user can set the 'budding_mother' parameter to control the fecundity of mother lineages. If 'budding_mother' is NA, then mother lineages will appear stochastic given the global probability of budding speciation events controlled by the 'budding_prob' argument. If 'budding_mother' is 0.0, then mother lineages will have only a single daughter lineage and the next speciation event will always be a symmetric speciation. If 'budding_mother' is set to 1.0, then once a budding speciation happens the surviving will always survive the speciation event and the descendant nodes will never represent symmetric speciation events. Intermediate values can be used to set intermediate scenarios.
#'
#' @param tree phylogeny with branch lengths.
#' @param Q transition matrix for a Markov model.
#' @param anc state at the root. If "NULL", then a random state is selected.
#' @param budding_prob probability that the speciation event was budding.
#' @param budding_mother autocorrelation of budding events. If NA, then no autocorrelation happens and the occurrence of budding speciation will always follow the same probability of 'budding_prob'. If 'budding_mother' is non-zero, then the probability that the mother lineage will generate another new lineage through a budding speciation event (and, thus, survive another speciation event) will be given by 'budding_mother' instead. The probability of the next budding event after a symmetrical speciation event will reset to 'budding_prob'.
#' @param change_rate the rate for an exponential reduction of the rate of trait evolution proportional to lineage-age.
#' @param decay_fn TRUE or FALSE. If TRUE, then the rate of trait evolution has a negative (exponential) relationship with lineage-age. If FALSE, then this relationship becomes positive.
#' @param cladogenetic_change if a trait change should occur in the node associated with a budding event. Option are "none": no change, anagenetic changes only; "flat": trait changes randomly; "prob": trait changes following the transition probabilities from the Q matrix.
#'
#' @return A list with simmap, tip_state, edge_state, ancestry, rate_scaler, and scaler_mat.
#' @export
#' @importFrom ape vcv.phylo reorder.phylo
#' @importFrom expm expm
sim_Mk_budding_exp <- function(tree, Q, anc = NULL, budding_prob = 0.0, budding_mother = NA, change_rate = 2.0, decay_fn = TRUE, cladogenetic_change = "none"){

    ## Check for the use of the budding_prob argument.
    if( budding_prob > 1.0 ) stop( "budding_prob need to be smaller than 1.0" )
    if( budding_prob < 0.0 ) stop( "budding_prob need to be larger than 0.0" )

    ## Check if budding_mother is set to NA.
    if( is.na( budding_mother ) ){
        budding_mother <- budding_prob
    } else{ ## Check if numeric.
        if( !is.numeric( budding_mother ) ) stop( "budding_mother needs to be NA or a value between 0.0 and 1.0" )
    }

    ## Check for the use of the budding_mother argument.
    if( budding_mother > 1.0 ) stop( "budding_mother need to be smaller than 1.0" )
    if( budding_mother < 0.0 ) stop( "budding_mother need to be larger than 0.0" )

    ## If no budding in the model, then budding_mother should also be absent.
    if( budding_prob == 0.0 ) budding_mother <- 0.0

    ## Check for the cladogenetic options:
    cladogenetic_change <- match.arg(arg = cladogenetic_change
                                     , choices=c("none","flat", "prob")
                                     , several.ok = FALSE)

    ## Define the size of the chunks to be used for the simulation.
    ## Max age helps to deal with non-ultrametric trees.
    chunk_length <- max( diag( vcv.phylo(tree) ) ) / 1000 ## 1000 pieces of tree age.

    ## Get the names for the states or return error message if fail.
    if( is.null(rownames(Q)) ){
        if( is.null(colnames(Q)) ){
            stop( "Q needs rownames or colnames equal to the states in the data." )
        } else{
            ss <- colnames(Q)
        }
    } else{
        ss <- rownames(Q)
    }

    ## Check if ancestral state is within states on the Q matrix:
    ## Set the root value for the simulation based on the arguments for the function.
    if( !is.null(anc) ){
        options_state <- ss
        a <- match.arg(arg = anc, choices = options_state, several.ok = FALSE)
    } else{
        a <- sample(ss, 1)
    }
    tt <- reorder.phylo(tree) ## reorder the tree cladewise.
    ## P <- vector(mode = "list", length = nrow(tt$edge))
    ## Need a 'maps' list in the same format as 'phytools' 'simmap' object.
    maps <- vector(mode = "list", length = nrow(tt$edge))

    ## For each edge of the tree we need an starting (corners) and an ending state (node).
    STATES <- matrix(NA, nrow(tt$edge), 2)

    ## For each edge of the tree we need a number to keep track of the lineage.
    ANCESTRY <- vector(mode = "numeric", length = nrow(tt$edge) )

    ## The rate scaler is a list in the same format of maps. So we can plot the rate scalers per branch.
    ## SCALER <- vector(mode = "numeric", length = nrow(tt$edge) )
    SCALER <- vector(mode = "list", length = nrow(tt$edge))

    root_id <- which(tt$edge[, 1] == Ntip(tt) + 1) ## The root edges.
    ## Set the root state for the simulation.
    STATES[root_id, 1] <- a
    ## Initiate the ancestry of the first lineage. At the edges connecting to the root we have the first and the second lineages
    ANCESTRY[root_id] <- c(1,2)
    ## The root have rate scaler of 1.0 because the history of the clade has just started. The rate scaler just after speciation (of a new lineage) will also be 1.0
    ## Note that 1.0 will be the highest scaler when decay_fn=TRUE and the lowest when decay_fn=FALSE. But in either case the starting scaler will be the same.
    SCALER[root_id] <- 1.0 ## These are the root lineages, they will start with 1.0 as the scaler.
    ## Starting point for the rate scaler is a constant.
    start_scaler <- 1.0

    ## Create a matrix for the storage of the rate scalers associated with the age.
    mat_scaler <- matrix(data = NA, nrow = 1, ncol = 4)

    for (j in 1:nrow(tt$edge)){

        ## For each edge we compute the exponential of the Q matrix multiplied by the time interval.
        ## These give the probabilities of transition after the branch length time has passed.
        ## When using an age_fn, the probability of transition at each time interval need to respond to the time.
        ## The time is the average of the lineage age between the start and the end of the interval (simulation chunk).

        ## ###################################
        ## SIMULATION ALONG THE BRANCH - ANAGENETIC CHANGES.
        ## This need to be done by the chunk.
        chunk_vec <- get_chunk_vec(ll = tt$edge.length[j], chunk_length = chunk_length)
        chunk_state <- vector(mode = "character", length = length(chunk_vec))
        chunk_scaler <- vector(mode = "numeric", length = length(chunk_vec))
        start_state <- STATES[j,1]

        ## Define the scaler, if we are in the root lineage or in a new lineage, then the start state of the scaler is 1.0
        ## This will need to change depending on the type of model. Each model will have a step function and a starting state.
        if( sum(ANCESTRY == ANCESTRY[j]) == 1 ){ ## First time the lineage number appears in the ANCESTRY vector.
            ## New lineage, including the root edges.
            new_lineage <- TRUE
        } else{
            ## Mother lineage. Lineage already exists.
            ## The ancestral edge NEEDS to belong to the same continued lineage.
            ancestral_edge <- tt$edge[,2] == tt$edge[j,1]
            new_lineage <- FALSE
        }

        for( w in 1:length(chunk_vec) ){
            ## Implement the lineage age function:
            chunk_sum <- ifelse(test = w == 1, yes = 0, no = sum(chunk_vec[1:(w-1)]) )
            ## get_chunk_age will compute the age of the current chunk following the trajectory of the current lineage.
            age_tmp <- get_chunk_age(tree = tt, current_node = tt$edge[j,1], chunk_length = chunk_vec[w]
                                     , chunk_sum = chunk_sum, lineage_number = ANCESTRY[j]
                                     , ancestry = ANCESTRY)
            ## Take one step on the scaler based on the change_rate and the value of the previous scaler.
            ## If the change_rate is 0.0, then this scaler will always be 1.0 and the model will be time-homogeneous.
            if( decay_fn ){
                age_tmp_scaler <- start_scaler / exp(change_rate * age_tmp)
            } else{
                ## Implements the positive relationship with time. Also exponential.
                age_tmp_scaler <- start_scaler * exp(change_rate * age_tmp)
            }
            ## Scale the Q matrix and draw the new state after the chunk time.
            Q_tmp <- expm(Q * chunk_vec[w] * age_tmp_scaler)
            rownames(Q_tmp) <- colnames(Q_tmp) <- ss
            new_state <- ss[ which( rmultinom(n = 1, size = 1, prob = Q_tmp[start_state,] )[, 1] == 1) ]
            ## Log the state and the scaler values.
            mat_scaler <- rbind(mat_scaler, c(new_lineage, age_tmp, chunk_vec[w], age_tmp_scaler))
            chunk_state[w] <- new_state
            start_state <- new_state
            chunk_scaler[w] <- age_tmp_scaler
        }

        ## Update the maps for the state and the scaler values.
        SCALER[[j]] <- data.frame(chunk = chunk_vec, scaler = chunk_scaler)
        maps[[j]] <- setNames(object = chunk_vec, nm = chunk_state) ## The maps list for the stochastic maps.
        new <- chunk_state[ length(chunk_state) ] ## The state at the end of the edge.
        ## ###################################

        ## ###################################
        ## STEPS AT THE DESCENDANT NODE - BUDDING SPECIATION.
        ## The state that was sampled now is the state at the descendant node.
        STATES[j, 2] <- new
        ## Check if the descendant node is an internal node or tip.
        ii <- which(tt$edge[, 1] == tt$edge[j, 2])
        if(length(ii) > 0){ ## Internal node (two descendant edges).

            ## #####################################
            ## Check for mother lineage and if budding happens.
            if( !new_lineage ){
                ## Apply 'budding_mother' parameter.
                budd_happens <- sample(x = c(TRUE, FALSE), size = 1, prob = c(budding_mother, 1-budding_mother) )
            } else{
                ## Now lineage, so apply 'budding_prob'.
                budd_happens <- sample(x = c(TRUE, FALSE), size = 1, prob = c(budding_prob, 1-budding_prob) )
            }
            ## #####################################

            ## #####################################
            ## Apply budding speciation if necessary.
            if( budd_happens ){
                ## One lineage keeps the same ancestry (at random).
                ndrd <- sample(x = c(1,2), size = 2, replace = FALSE)
                ## Update the continuation of the ancestral.
                ANCESTRY[ii[ndrd[1]]] <- ANCESTRY[j] ## j is the ancestral edge.
                STATES[ii[ndrd[1]], 1] <- new
                ## Update the new lineage, with the cladogenetic function.
                ANCESTRY[ii[ndrd[2]]] <- max(ANCESTRY) + 1 ## the new lineage.

                ## Check if we make a cladogenetic change to the trait:
                if( cladogenetic_change == "flat" ){
                    ## Have a deterministic change in the state at the budding speciation event.
                    ## If binary trait, then bounce to the other state.
                    budd_states <- ss[ ss != new ] ## Any of the other states.
                    STATES[ii[ndrd[2]], 1] <- sample(x = budd_states, size = 1)
                } else if( cladogenetic_change == "prob" ){
                    ## Change the trait state using the prob from the Q transition matrix.
                    row_sub <- ss == new
                    col_sub <- !row_sub
                    prob_shift <- Q[row_sub, col_sub] / sum( Q[row_sub, col_sub] )
                    STATES[ii[ndrd[2]], 1] <- sample(x = ss[col_sub], prob = prob_shift, size = 1)
                } else{
                    ## Don't change the trait.
                    STATES[ii[ndrd[2]], 1] <- new
                }

            } else{
                ## Both descendants are new lineages and they share the same state.
                STATES[ii, 1] <- new
                ANCESTRY[ii] <- c(1,2) + max( ANCESTRY )
            }

        }
        ## ###################################
    }

    xx_temp <- sapply(1:Ntip(tt), function(n, S, E) S[which(E == n)], S = STATES[, 2], E = tt$edge[,2])
    x <- as.factor( setNames(object = xx_temp, nm = tt$tip.label) )

    ## Make the mapped.edge to complete to simmap object.
    for(w in 1:length(maps) ){
        for(k in 1:length(maps[[w]]) ){
            mapped.edge[w,names(maps[[w]])[k]] <- mapped.edge[w,names(maps[[w]])[k]] + maps[[w]][k]
        }
    }

    ## Include all the simmap information in the phylogeny to return.
    tt$maps <- maps
    tt$Q <- Q
    tt$mapped.edge <- mapped.edge
    tt$logL <- 42 ## Any number would do here.
    class( tt ) <- c( class( tt ), "simmap" )

    ## Temporary object to return the rate scalers.
    colnames( mat_scaler ) <- c("new_lineage", "age", "chunk_length", "scaler")

    return( list(simmap = tt, tip_state = x, edge_state = STATES
                 , ancestry = ANCESTRY, rate_scaler_maps = SCALER
                 , scaler_mat = mat_scaler) )
}

#' Extract budding history from simulations
#'
#' This function can be used to extract the history of budding events across the branches and nodes of the tree.
#'
#' @param sim A simulation object generated with 'sim_Mk_budding_exp', 'sim_Mk_trace_history', 'sim_BM_budding_exp', or 'sim_BM_trace_history'.
#'
#' @return a list. budd_nodes are the nodes representing budding speciation events. budd_edges is an index, in the same order of phy$edge, selecting all branches belonging to mother lineages. A mother lineage is a lineage that survived one or more speciation events (because of budding).
#' @export
get_Budding_History <- function(sim){
    ## The budd_nodes show each node that a budding speciation event happened.
    ## The budd_edges marks all edges with a mother lineage (lineages that survived at least one speciation event.
    budd_lineages <- unique( sim$ancestry[duplicated(sim$ancestry)] )
    budd_nodes <- vector(mode = "numeric")
    budd_edges <- rep(x = FALSE, times = length(sim$ancestry) )
    if( "simmap" %in% names(sim) ){
        edge_mat <- sim$simmap$edge
    } else if( "contsim" %in% names(sim) ){
        edge_mat <- sim$contsim$edge
    } else{
        stop( "Check if 'sim' is the correct object!" )
    }
    for( i in budd_lineages ){
        anc_id <- sim$ancestry == i
        budd_edges <- anc_id | budd_edges ## This will update the budd_edges correctly.
        budd_id_nd <- duplicated( c( edge_mat[anc_id,] ) )
        budd_nodes <- append(x = budd_nodes
                             , values = c( edge_mat[anc_id,] )[budd_id_nd]
                             , after = length(budd_nodes) )
    }
    return( list(budd_nodes = budd_nodes, budd_edges = budd_edges) )
}

## A simple function to create a vector with the chunks.
##' @noRd
get_chunk_vec <- function(ll, chunk_length){
    chunk_div <- ll / chunk_length
    chunk_remainder <- chunk_div - floor(x = chunk_div)
    ## The last chunk will be just a proportion of the chunk_size.
    last_chunk <- chunk_remainder * chunk_length
    chunk_n <- ceiling(chunk_div)
    chunk_vec <- rep(x = chunk_length, length = chunk_n)
    chunk_vec[chunk_n] <- last_chunk
    return( chunk_vec )
}

##' @noRd
##' @importFrom ape node.depth.edgelength
get_chunk_age <- function(tree, current_node, chunk_length, chunk_sum, lineage_number, ancestry){
    ## Compute the age of the lineage so far.
    ## This is in the same unit as the edge lengths of the tree.
    tree_age_nodes <- round(node.depth.edgelength(phy = tree)[1], digits = 7) - round(node.depth.edgelength(phy = tree), digits = 7)
    lineage_nodes <- tree$edge[which(ancestry == lineage_number), 1] ## One of these, the oldest.
    age_ancestral <- max( tree_age_nodes[ lineage_nodes ] )
    age_current_node <- tree_age_nodes[ current_node ]
    ## The time interval of the lineage plus the length of the branch visited so far.
    age_lineage <- (age_ancestral - age_current_node) + chunk_sum + (chunk_length/2)
    return( age_lineage )
}

## Get lineages longer than 1 branch and create palette of colors.
## This vector of colors can be used to highlight the mother lineages.
##' @noRd
get_edge_color_lineages <- function(simMK, base_color = "gray"){
    mother_lineages <- as.numeric( which( table( simMK$ancestry ) > 1 ) )
    color_pool <- colors(distinct = TRUE)
    ## Get only basic colors, strip the numbers.
    color_pool <- unique( gsub(pattern = "[[:digit:]]+", replacement = "", x = color_pool) )
    ## Exclude the base color and any of its variants.
    color_pool <- color_pool[ -grep(pattern = base_color, x = color_pool) ]
    if( length(mother_lineages) > length( color_pool ) ){
        ## Protect for really large phylogenies.
        mult_color <- ceiling( length(mother_lineages) / length( color_pool ) )
        mother_palette <- c( sapply(1:mult_color, function(i) sample(x = color_pool, size = length(color_pool), replace = FALSE) ) )
    } else{
        mother_palette <- sample(x = color_pool, size = length(mother_lineages), replace = FALSE)
    }
    edge_colors <- rep(base_color, times = nrow(simMK$simmap$edge) )
    for( mm in 1:length(mother_lineages) ){
        edge_colors[ simMK$ancestry == mother_lineages[mm] ] <- mother_palette[mm]
    }
    return( edge_colors )
}

##' @noRd
get_des_scaler <- function(scaler, id){
    ## Get the last rate scaler value in the ancestral lineage.
    num_in <- which( id )
    x_mat <- scaler[[num_in]]
    last <- nrow( x_mat )
    x_last <- as.numeric( x_mat[last,2] )
    return( x_last )
}

## Function to plot the phylogeny and mark the mother lineages.
## Mother lineages are the lineages that survived through one or more events of speciation.
#' Plot phylogeny and show mother lineages
#'
#' Function plots the phylogenetic tree from a buddPhy simulation and highlights the mother lineages. Mother lineages are those that survived through one of more events of speciation via budding.
#' At the moment the function use colors randomly and, unfortunately, the color scheme cannot be chosen. Please re-run the plot to change color configuration.
#' @param sim_obj the output of a buddPhy simulation function.
#' @param background_color the color of the non-mother lineages.
#' @param edge_width the width of the phylogeny branches.
#' @param no_margin if to add margin. See 'ape::plot.phylo'.
#'
#' @return Plots a phylogenetic tree.
#' @export
#' @importFrom ape plot.phylo
plot_mother_lineages <- function(sim_obj, background_color = "gray", edge_width = 3, no_margin = TRUE){
    edge_colors <- get_edge_color_lineages(simMK = sim_obj, base_color = background_color)
    plot.phylo(x = sim_obj$simmap, edge.color = edge_colors, edge.width = edge_width, no.margin = no_margin)
}

#' Make stochastic map with rate scalers
#'
#' Plot the rate scalers used in the simulation to the phylogenetic tree. Also return a vector of colors to be used for visualization of the gradient paired with the breaking points for the rate scaler values. \cr
#' Use 'low_cut' or 'high_cut' to focus the rate categories into the slow or fast rate categories. This can be used to provide a more informative plot of rates along the branches of the phylogeny.
#'
#' @param sims the simulation object output from 'sim_Mk_budding_exp'.
#' @param ncat number of categories to make the legend for the scalers plot.
#' @param low_cut break for the first category (slow rate scalers). Rate scaler values above 'low_cut' value will be divided into 'ncat'-1 categories.
#' @param high_cut break for the last category (fast rate scalers). Rate scaler values below 'high_cut' value will be divided into 'ncat'-1 categories.
#'
#' @return a list with a simmap object and a matrix with legend information
#' @export
make_scaler_simmap <- function(sims, ncat = 10, low_cut, high_cut){

    scaler_maps <- sims$rate_scaler_maps
    raw_poll <- unlist( sapply(scaler_maps, function(x) x[,2]) )
    scaler_pool <- range( raw_poll )

    if( missing( high_cut ) & missing( low_cut ) ){
        ## Compute natural breaks.
        breaks_right <- seq(from = scaler_pool[1], to = scaler_pool[2], length.out = ncat+1)
        breaks_right[1] <- floor( breaks_right[1] )
        breaks_right[length(breaks_right)] <- ceiling( breaks_right[length(breaks_right)] )
    }

    if( !missing( high_cut ) & !missing( low_cut ) ) stop("can only provide low_cut OR high_cut, not both.")

    if( !missing( high_cut ) ){
        ## Use and check the value of high_cut.
        if( !is.numeric( high_cut ) ) stop("if provided, high_cut needs to be a numeric value.")
        if( high_cut > max(scaler_pool) ) stop("high_cut value is larger than the max observed scaler.")
        breaks_right <- seq(from = scaler_pool[1], to = high_cut, length.out = ncat)
        ## Need to fix both ends to make sure the breaks include the observed values.
        breaks_right <- c(breaks_right, ceiling(scaler_pool[2]))
        breaks_right[1] <- floor( breaks_right[1] )
    }

    if( !missing( low_cut ) ){
        ## Use and check the value of low_cut.
        if( !is.numeric( low_cut ) ) stop("if provided, low_cut needs to be a numeric value.")
        if( low_cut < min(scaler_pool) ) stop("low_cut value is smaller than the min observed scaler.")
        breaks_right <- seq(from = low_cut, to = scaler_pool[2], length.out = ncat)
        ## Need to fix both ends to make sure the breaks include the observed values.
        breaks_right <- c(floor(scaler_pool[1]), breaks_right)
        breaks_right[length(breaks_right)] <- ceiling( breaks_right[length(breaks_right)] )
    }

    cat_labels <- LETTERS[1:ncat]
    maps_formated <- vector(mode = "list", length = length(scaler_maps))
    for( i in 1:length(scaler_maps) ){
        maps_labels <- as.character( cut(x = scaler_maps[[i]][,2], breaks = breaks_right
                                         , labels = cat_labels, right = TRUE) )
        maps_formated[[i]] <- setNames(object = scaler_maps[[i]][,1], nm = maps_labels)
    }

    ## Make the mapped.edge to complete to simmap object.
    mapped.edge <- matrix(0, nrow = nrow(sims$simmap$edge), ncol = ncat, dimnames = list(paste(sims$simmap$edge[,1],",",sims$simmap$edge[,2],sep=""), cat_labels))
    for(w in 1:length(maps_formated) ){
        for(k in 1:length(maps_formated[[w]]) ){
            mapped.edge[w,names(maps_formated[[w]])[k]] <- mapped.edge[w,names(maps_formated[[w]])[k]] + maps_formated[[w]][k]
        }
    }

    ## Make a color vector for plot simmap.
    colorfn <- colorRampPalette(colors = c("green", "yellow", "red"))
    col_rgb <- colorfn(ncat)
    col_vec <- setNames(object = col_rgb, nm = cat_labels)

    ## Create a matrix with the legend for the plot.
    legend_matrix <- data.frame(from = breaks_right[-length(breaks_right)], to = breaks_right[-1], label = cat_labels, col = col_rgb)

    ## Substitute the maps from the original tree:
    simmap <- sims$simmap
    simmap$maps <- maps_formated
    simmap$mapped.edge <- mapped.edge

    ## Return the stochastic maps and the matrix with the legend.
    ll <- list(simmap = simmap, legend_matrix = legend_matrix, col_vec = col_vec)
    return( ll )
}
