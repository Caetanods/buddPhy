#' Simulate discrete trait evolution tracing lineage history.
#'
#' This function simulates discrete traits based on the same lineage history created when using 'sim_Mk_budding_exp'. It is useful to simulate under the same lineage history and phylogenetic tree while changing rate parameters or parameter for the decay function. If you don't want to retain a particular lineage history, then please refer to 'sim_Mk_budding_exp'.
#'
#' @param sim_MK the output from ''sim_Mk_budding_exp'.
#' @param Q transition matrix for a Markov model.
#' @param anc state at the root. If "NULL", then a random state is selected.
#' @param change_rate the rate for an exponential reduction of the rate of trait evolution proportional to lineage-age.
#' @param decay_fn TRUE or FALSE. If TRUE, then the rate of trait evolution has a negative (exponential) relationship with lineage-age. If FALSE, then this relationship becomes positive.
#' @param cladogenetic_change if a trait change should occur in the node associated with a budding event. Option are "none": no change, anagenetic changes only; "flat": trait changes randomly; "prob": trait changes following the transition probabilities from the Q matrix.
#'
#' @return A list with simmap, tip_state, edge_state, ancestry, rate_scaler, and scaler_mat.
#' @export
#' @importFrom ape vcv.phylo reorder.phylo
#' @importFrom expm expm
#' @importFrom ratematrix mergeSimmap
sim_Mk_trace_history <- function(sim_MK, Q, anc = NULL, change_rate = 2.0, decay_fn = TRUE, cladogenetic_change = "none"){

    ## Check for the cladogenetic options:
    cladogenetic_change <- match.arg(arg = cladogenetic_change
                                     , choices=c("none","flat", "prob")
                                     , several.ok = FALSE)

    ## Define the size of the chunks to be used for the simulation.
    ## Max age helps to deal with non-ultrametric trees.
    tt <- mergeSimmap(phy = sim_MK$simmap, drop.regimes = TRUE)
    chunk_length <- max( diag( vcv.phylo(tt) ) ) / 1000 ## 1000 pieces of tree age.

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
    ## Need a 'maps' list in the same format as 'phytools' 'simmap' object.
    maps <- vector(mode = "list", length = nrow(tt$edge))

    ## For each edge of the tree we need an starting (corners) and an ending state (node).
    STATES <- matrix(NA, nrow(tt$edge), 2)

    ## For each edge of the tree we need a number to keep track of the lineage.
    ANCESTRY <- vector(mode = "numeric", length = nrow(tt$edge) )

    ## However, here we already know the history of all lineages. We need to track this history instead of creating a new one. So here we define DESTINY:
    DESTINY <- sim_MK$ancestry ## The object from the previous simulation.

    ## Create a vector with the complete history of budding per edge.
    BUDD_HISTORY <- vector(mode = "logical", length = nrow(tt$edge))
    for( i in 1:nrow(tt$edge) ){
        ## Need to track the descendants. If any of the descendants is a continuation of the current lineage, then we need to set this edge as a budding edge.
        des_node <- tt$edge[i,2] ## Want to know if at end of this edge we have a budding event.
        desc_lineages <- DESTINY[ tt$edge[,1] == des_node ]
        curr_lineage <- DESTINY[i]
        BUDD_HISTORY[i] <- curr_lineage %in% desc_lineages
    }

    ## The rate scaler is a list in the same format of maps. So we can plot the rate scalers per branch.
    SCALER <- vector(mode = "list", length = nrow(tt$edge))

    root_id <- which(tt$edge[, 1] == Ntip(tt) + 1) ## The root edges.
    ## We always refer to the DESTINY object.
    ANCESTRY[root_id] <- DESTINY[root_id]
    ## Set the root state for the simulation.
    STATES[root_id, 1] <- a
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
            ## start_scaler <- 1.0
            new_lineage <- TRUE
        } else{
            ## Mother lineage. Lineage already exists.
            ## The ancestral edge NEEDS to belong to the same continued lineage.
            ancestral_edge <- tt$edge[,2] == tt$edge[j,1]
            ## start_scaler <- get_des_scaler(scaler = SCALER, id = ancestral_edge)
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

            ## Need to know if a budding event happened at this node.
            budd_happens <- BUDD_HISTORY[j]

            ## #####################################
            ## Apply budding speciation if necessary.
            if( budd_happens ){
                ## Find which lineage keeps the same ancestry.
                mother_lineage_id <- ii[ DESTINY[j] == DESTINY[ ii ] ]
                daughter_lineage_id <- ii[ DESTINY[j] != DESTINY[ ii ] ]
                ## Updates the ancestry, here with its destiny.
                ANCESTRY[ii] <- DESTINY[ii]
                ## The mother lineage keeps the same state:
                STATES[mother_lineage_id, 1] <- new

                ## Check if we make a cladogenetic change to the trait:
                if( cladogenetic_change == "flat" ){
                    ## Have a deterministic change in the state at the budding speciation event.
                    ## If binary trait, then bounce to the other state.
                    budd_states <- ss[ ss != new ] ## Any of the other states.
                    STATES[daughter_lineage_id, 1] <- sample(x = budd_states, size = 1)
                } else if( cladogenetic_change == "prob" ){
                    ## Change the trait state using the prob from the Q transition matrix.
                    row_sub <- ss == new
                    col_sub <- !row_sub
                    prob_shift <- Q[row_sub, col_sub] / sum( Q[row_sub, col_sub] )
                    STATES[daughter_lineage_id, 1] <- sample(x = ss[col_sub], prob = prob_shift, size = 1)
                } else{
                    ## Don't change the trait. No special event at cladogenesis.
                    STATES[daughter_lineage_id, 1] <- new
                }

            } else{
                ## Both descendants are new lineages and they share the same state.
                STATES[ii, 1] <- new
                ANCESTRY[ii] <- DESTINY[ii]
            }

        }
        ## ###################################
    }

    xx_temp <- sapply(1:Ntip(tt), function(n, S, E) S[which(E == n)], S = STATES[, 2], E = tt$edge[,2])
    x <- as.factor( setNames(object = xx_temp, nm = tt$tip.label) )

    ## Make the mapped.edge to complete to simmap object.
    mapped.edge <- matrix(0, nrow = nrow(tt$edge), ncol = ncol(Q)
                          , dimnames = list(paste(tt$edge[,1],",",tt$edge[,2],sep=""), ss))
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
