## Many of these functions are based on Liam's code for stochastic mapping in phytools.

## How should we compute the chunk length? Chunk length is only necessary for models with some age_fn. If Q is constant, then we can use the exp to get the waiting time until the next event (just like a standard MK model). If Q is going to change at each chunk following some age_fn, then we need to simulate the trait at time t after delta_t. Meaning that each chunk will be the branch between two pseudo-nodes and will have some rate multiplier applied to that chunk (controlled by the age_fn).

## #######################################################
## EXPANDING THE FUNCTION
## The type of decay function and/or increase in rate of trait evolution with lineage-age can be expanded. Here I will create one function for each type of lineage-age relationship. It is simpler.
## Need to add a probability of autocorrelation among the budding events. How likely is a budding event to follow another one.

## Function to simulate discrete traits under budding speciation.
## Multiple factors control the simulation and the probability of budding.
## Function keeps track of the ancestral and descendant lineages.
## Ideally we want a general function that can take any age function to make the simulation.
## Here we implemented an exponential decay function with the "change_rate" parameter. This simulates an exponential negative relationship between lineage-age and the rate of trait evolution.

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
sim_Mk_budding_exp <- function(tree, Q, anc = NULL, budding_prob = 0.0, budding_mother = NA, change_rate = 2.0, decay_fn = TRUE, cladogenetic_change = "none"){

    ## Next thing to do: Add an autocorrelation value for the budding process. [working on this]
    ## Implement other lineage_age processes (in this or other functions).

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
    cladogenetic_change <- match.arg(arg = cladogenetic_change, choices=c("none","flat", "prob"), several.ok = FALSE)

    ## Check if ancestral state is within states on the Q matrix:
    if( !is.null(anc) ){
        options_state <- rownames(Q)
        anc <- match.arg(arg = anc, choices = options_state, several.ok = FALSE)
    }

    ## Define the size of the chunks to be used for the simulation.
    chunk_length <- vcv.phylo(tree)[1,1] / 1000 ## 1000 pieces of tree age.

    ss <- rownames(Q)
    tt <- reorder.phylo(tree) ## reorder the tree cladewise.
    P <- vector(mode = "list", length = nrow(tt$edge))
    ## Need a 'maps' list in the same format as 'phytools' 'simmap' object.
    maps <- vector(mode = "list", length = nrow(tt$edge))

    ## Set the root value for the simulation based on the arguments for the function.
    if (is.null(anc)){
        a <- sample(ss, 1)
    } else{
        a <- anc
    }

    ## For each edge of the tree we need an starting and an ending state.
    STATES <- matrix(NA, nrow(tt$edge), 2)
    ## For each edge of the tree we need a number to keep track of the lineage.
    ANCESTRY <- vector(mode = "numeric", length = nrow(tt$edge) )
    ## Need to keep track of the rate scaler to apply the exponential decay of the rate.
    ## The rate scaler is populated at the end of the simulation on the branch. The descendant branch might refer to this value if the speciation was a budding speciation event. Because the age of the lineage accumulates.
    SCALER <- vector(mode = "numeric", length = nrow(tt$edge) )
    root_id <- which(tt$edge[, 1] == Ntip(tt) + 1) ## The root edges.
    ## Set the root state for the simulation.
    STATES[root_id, 1] <- a
    ## Initiate the ancestry of the first lineage. At the edges connecting to the root we have the first and the second lineages
    ANCESTRY[root_id] <- c(1,2)
    ## The root have rate scaler of 1.0 because the history of the clade has just started. The rate scaler just after speciation (the new lineage) will also be 1.0
    SCALER[root_id] <- 1.0 ## These are the root lineages, they will start with 1.0 as the scaler.

    ## Create a matrix for the storage of the rate scalers associated with the age.
    mat_scaler <- matrix(data = NA, nrow = 1, ncol = 3)

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
        start_state <- STATES[j,1]

        ## Define the scaler, if we are in the root lineage or in a new lineage, then the start state of the scaler is 1.0
        ## This will need to change depending on the type of model. Each model will have a step function and a starting state.
        if( sum(ANCESTRY == ANCESTRY[j]) == 1 ){ ## First time the lineage number appears in the ANCESTRY vector.
            ## New lineage, including the root edges.
            start_scaler <- 1.0
            new_lineage <- TRUE
        } else{
            ## Mother lineage. Lineage already exists.
            ## The ancestral edge NEEDS to belong to the same continued lineage.
            ancestral_edge <- tt$edge[,2] == tt$edge[j,1]
            start_scaler <- SCALER[ancestral_edge]
            new_lineage <- FALSE
        }

        for( w in 1:length(chunk_vec) ){
            ## Implement the lineage age function:
            chunk_sum <- ifelse(test = w == 1, yes = 0, no = sum(chunk_vec[1:(w-1)]) )
            ## get_chunk_age will compute the age of the current chunk following the trajectory of the current lineage.
            age_tmp <- get_chunk_age(tree = tt, current_node = tt$edge[j,1], chunk_length = chunk_vec[w]
                                         , chunk_sum = chunk_sum, lineage_number = ANCESTRY[j], ancestry = ANCESTRY)
            ## Take one step on the scaler based on the change_rate and the value of the previous scaler.
            ## If the change_rate is 0.0, then this scaler will always be 1.0 and the model will be time-homogeneous.
            if( decay_fn ){
                start_scaler <- start_scaler / exp(change_rate * age_tmp)
            } else{
                ## Implements the positive relationship with time. Also exponential.
                start_scaler <- start_scaler * exp(change_rate * age_tmp)
            }
            ## Keep a matrix with all the scalers used in the simulation. [This is just to follow the flow of the simulations.]
            mat_scaler <- rbind(mat_scaler, c(new_lineage, age_tmp, start_scaler))
            ## The scaler is lower than 1.0, so the rate will decrease with the age of the lineage.
            Q_tmp <- expm(Q * chunk_vec[w] * start_scaler)
            rownames(Q_tmp) <- colnames(Q_tmp) <- ss
            new_state <- ss[ which( rmultinom(n = 1, size = 1, prob = Q_tmp[start_state,] )[, 1] == 1) ]
            chunk_state[w] <- new_state
            start_state <- new_state
        }

        ## Update the rate scaler at the end of the branch for this branch.
        SCALER[j] <- start_scaler
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
                    row_sub <- rownames(Q) == new
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
    mapped.edge <- matrix(0, nrow = nrow(tt$edge), ncol = ncol(Q), dimnames = list(paste(tt$edge[,1],",",tt$edge[,2],sep=""), rownames(Q)))
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
    colnames( mat_scaler ) <- c("new_lineage", "age", "scaler")

    return( list(simmap = tt, tip_state = x, edge_state = STATES
                 , ancestry = ANCESTRY, rate_scaler = SCALER
                 , scaler_mat = mat_scaler) )
}

## Function matches the edge matrix of the phylogeny with the states matrix from the budding simulation.
## Then returns the nodes where budding speciation was simulated to have happened.
##' @noRd
get_budding_nodes <- function(simMK){
    tree <- simMK$simmap
    phy_edge <- tree$edge
    state_edge <- simMK$edge_state
    root_node <- Ntip( tree ) + 1
    max_node <- max( phy_edge[,1] )
    budding_node <- vector(mode = "logical", length = length(root_node:max_node) )
    nodes <- root_node:max_node
    for( nd in 1:length(nodes) ){
        node_pos <- which( phy_edge[,1] == nodes[nd] )
        ## If a single occurrence, then 0. If 2 or more occurrences (therefore the lineage survived speciation) then it is >0 (transformed into TRUE).
        budding_node[nd] <- as.logical( length( unique( state_edge[node_pos,1] ) ) - 1 )
    }
    return( setNames(object = budding_node, nm = nodes) )
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
    mother_palette <- sample(x = color_pool, size = length(mother_lineages), replace = FALSE)
    edge_colors <- rep(base_color, times = nrow(simMK$simmap$edge) )
    for( mm in 1:length(mother_lineages) ){
        edge_colors[ simMK$ancestry == mother_lineages[mm] ] <- mother_palette[mm]
    }
    return( edge_colors )
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
    bud_nodes <- get_budding_nodes(simMK = sim_obj)
    edge_colors <- get_edge_color_lineages(simMK = sim_obj, base_color = background_color)
    plot.phylo(x = sim_obj$simmap, edge.color = edge_colors, edge.width = edge_width, no.margin = no_margin)
}
