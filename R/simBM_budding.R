#' Simulate continuous trait evolution with budding and exponential decay function.
#'
#' Similar to 'sim_Mk_budding_exp' function, but simulate continuous traits using a Brownian motion model instead of discrete traits.
#' \cr
#' This function adds a probability of budding speciation and allows only for a negative relationship between rates of trait evolution and lineage-age. Function also allows for cladogenetic change in the trait (a deterministic change in the state of the trait after budding speciation happens).
#' \cr
#' MOTHER LINEAGES: The user can set the 'budding_mother' parameter to control the fecundity of mother lineages. If 'budding_mother' is NA, then mother lineages will appear stochastic given the global probability of budding speciation events controlled by the 'budding_prob' argument. If 'budding_mother' is 0.0, then mother lineages will have only a single daughter lineage and the next speciation event will always be a symmetric speciation. If 'budding_mother' is set to 1.0, then once a budding speciation happens the surviving will always survive the speciation event and the descendant nodes will never represent symmetric speciation events. Intermediate values can be used to set intermediate scenarios.
#'
#' @param tree phylogeny with branch lengths.
#' @param sigma the base rate for the Brownian-motion model.
#' @param anc the ancestral state for the trait.
#' @param budding_prob probability that the speciation event was budding.
#' @param budding_mother autocorrelation of budding events. If NA, then no autocorrelation happens and the occurrence of budding speciation will always follow the same probability of 'budding_prob'. If 'budding_mother' is non-zero, then the probability that the mother lineage will generate another new lineage through a budding speciation event (and, thus, survive another speciation event) will be given by 'budding_mother' instead. The probability of the next budding event after a symmetrical speciation event will reset to 'budding_prob'.
#' @param change_rate the rate for an exponential reduction of the rate of trait evolution proportional to lineage-age.
#' @param decay_fn TRUE or FALSE. If TRUE, then the rate of trait evolution has a negative (exponential) relationship with lineage-age. If FALSE, then this relationship becomes positive.
#' @param jump_size controls the size of the cladogenetic change associated with budding speciation. The new lineage will have an average trait value shifted by a proportion equal to 'jump_size'. For example, 'jump_size = 0.5' will cause the new lineages to start with an average trait value equal to 50\% above or below the current trait value. Setting to NA (default) will turn-off the option of cladogenetic trait changes.
#'
#' @return A list with the information from the simulation and the states at the tips.
#' @export
#' @importFrom ape vcv.phylo reorder.phylo
sim_BM_budding_exp <- function(tree, sigma, anc = 0.0, budding_prob = 0.0, budding_mother = NA, change_rate = 2.0, decay_fn = TRUE, jump_size = NA){

    ## This function is very similar to the case with a discrete trait evolution model.
    ## We can either track the value of the trait at each branch length interval or we can simply rescale the tree following the rate function and perform a BM simulation afterwards.

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

    ## Define the size of the chunks to be used for the simulation.
    ## Max age helps to deal with non-ultrametric trees.
    chunk_length <- max( diag( vcv.phylo(tree) ) ) / 1000 ## 1000 pieces of tree age.

    ## Check if the jump_size is valid.
    if( is.numeric(jump_size) & jump_size < 0.0 ) stop("jump_size cannot be a negative number.")

    tt <- reorder.phylo(tree) ## reorder the tree cladewise.
    ## Need a 'maps' list in the same format as 'phytools' 'simmap' object.
    maps <- vector(mode = "list", length = nrow(tt$edge))

    ## For each edge of the tree we need an starting and an ending state.
    STATES <- matrix(NA, nrow(tt$edge), 2)
    ## For each edge of the tree we need a number to keep track of the lineage.
    ANCESTRY <- vector(mode = "numeric", length = nrow(tt$edge) )

    ## The rate scaler is a list in the same format of maps. So we can plot the rate scalers per branch.
    SCALER <- vector(mode = "list", length = nrow(tt$edge))

    root_id <- which(tt$edge[, 1] == Ntip(tt) + 1) ## The root edges.
    ## Set the root state for the simulation.
    STATES[root_id, 1] <- as.numeric( anc )
    ## Initiate the ancestry of the first lineage.
    ANCESTRY[root_id] <- c(1,2)
    SCALER[root_id] <- 1.0
    start_scaler <- 1.0

    ## Create a matrix for the storage of the rate scalers associated with the age.
    mat_scaler <- matrix(data = NA, nrow = 1, ncol = 4)

    for (j in 1:nrow(tt$edge)){

        ## ###################################
        ## SIMULATION ALONG THE BRANCH - ANAGENETIC CHANGES.
        ## This need to be done by the chunk.
        chunk_vec <- get_chunk_vec(ll = tt$edge.length[j], chunk_length = chunk_length)
        chunk_state <- vector(mode = "character", length = length(chunk_vec))
        chunk_scaler <- vector(mode = "numeric", length = length(chunk_vec))
        start_state <- as.numeric(STATES[j,1]) ## Continuous value.

        if( sum(ANCESTRY == ANCESTRY[j]) == 1 ){
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
            ## Scale the BM rate and simulate the trait evolution at the chunk.
            sigma_tmp <- sigma * chunk_vec[w] * age_tmp_scaler
            new_state <- as.numeric(start_state) + rnorm(n = 1, sd = sqrt(sigma_tmp) )
            ## Log the state and the scaler values.
            mat_scaler <- rbind(mat_scaler, c(new_lineage, age_tmp, chunk_vec[w], age_tmp_scaler))
            chunk_state[w] <- new_state
            start_state <- new_state ## Update the state at the chunk.
            chunk_scaler[w] <- age_tmp_scaler
        }

        ## Update the maps for the state and the scaler values.
        ## Note that the state here is a numeric value.
        SCALER[[j]] <- data.frame(chunk = chunk_vec, scaler = chunk_scaler)
        maps[[j]] <- setNames(object = chunk_vec, nm = chunk_state)
        new <- as.numeric( chunk_state[ length(chunk_state) ] ) ## The state at the end of the edge.
        ## ###################################

        ## ###################################
        ## STEPS AT THE DESCENDANT NODE - BUDDING SPECIATION.
        ## The state that was sampled now is the state at the descendant node.
        STATES[j, 2] <- as.numeric( new )
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

                ## Check if we should make a jump at the new species.
                if( is.na(jump_size) ){
                    STATES[ii[ndrd[2]], 1] <- new
                } else{
                    ## Make the cladogenetic jump.
                    STATES[ii[ndrd[2]], 1] <- new + (sample(x = c(1.0,-1.0), size = 1) * new * jump_size)
                }
            } else{
                ## Both descendants are new lineages and they share the same state.
                STATES[ii, 1] <- new
                ANCESTRY[ii] <- c(1,2) + max( ANCESTRY )
            }

        }
        ## ###################################
    }

    ## Get the states for the tips:
    xx_temp <- sapply(1:Ntip(tt), function(n, S, E) S[which(E == n)], S = STATES[, 2], E = tt$edge[,2])
    x <- setNames(object = as.numeric(xx_temp), nm = tt$tip.label)

    ## Temporary object to return the rate scalers.
    colnames( mat_scaler ) <- c("new_lineage", "age", "chunk_length", "scaler")

    return( list(contsim = tt, tip_state = x, edge_state = STATES
                 , ancestry = ANCESTRY, rate_scaler_maps = SCALER
                 , scaler_mat = mat_scaler, maps = maps) )
}

#' Simulate continuous trait evolution tracing lineage history.
#'
#' This function simulates continuous traits based on the same lineage history created when using 'sim_BM_budding_exp'. It is useful to simulate under the same lineage history and phylogenetic tree while changing rate parameters or parameter for the decay function. If you don't want to retain a particular lineage history, then please refer to 'sim_BM_budding_exp'.
#'
#' If 'jump_type' is set to "fixed", then each jump of character value at the node will be an increase or decrease of a fixed proportion of the value at that node. The proportion is controlled with the argument "jump_size". For example, 'jump_size = 0.5' will cause the new lineages to start with an average trait value equal to 50\% above or below the current trait value. On the other hand, if 'jump_type' is of type "gaussian", the jump in character state will be sampled from a normal distribution with mean equal to the current trait value and standard deviation (sd) equal to 'jump_size'.
#'
#' The parameter 'jump_all' controls the concept of new lineage. If set to TRUE, then a branching event which produces a new lineage will be followed by a jump on the phenotype value. For example, if a symmetric speciation event takes place, then both descendants are new lineages and both jump in the morphospace (controlled by the 'jump_size' and 'jump_type' parameters). If set to FALSE, then jumps in the morphospace associated with the origination of a new lineage occurs only if a budding speciation takes place.
#'
#' @param sim_BM output from the 'sim_BM_budding_exp'.
#' @param sigma the base rate for the Brownian-motion model.
#' @param anc the ancestral state for the trait.
#' @param change_rate the rate for an exponential reduction of the rate of trait evolution proportional to lineage-age.
#' @param decay_fn TRUE or FALSE. If TRUE, then the rate of trait evolution has a negative (exponential) relationship with lineage-age. If FALSE, then this relationship becomes positive.
#' @param jump_size controls the size of the cladogenetic change associated with budding speciation. Setting to NA (default) will turn-off the option of cladogenetic trait changes. If 'jump_type' is "fixed" then 'jump_size' is a proportion of the current trait value. If 'jump_type' is "gaussian", then 'jump_size' is the standard variation of a normal (Gaussian) distribution with mean equal to the current trait value and standard deviation equal to 'jump_size'.
#' @param jump_type if the phenotype jump at nodes should be "fixed" or "gaussian". See Details.
#' @param jump_all if jumps in the phenotype value should occur in all new lineages or only in new species originated from a mother lineage. See Details.
#'
#' @return A list with the information from the simulation and the states at the tips.
#' @export
#' @importFrom ape vcv.phylo reorder.phylo
sim_BM_trace_history <- function(sim_BM, sigma, anc = 0.0, change_rate = 2.0, decay_fn = TRUE, jump_size = NA, jump_type = "fixed", jump_all = FALSE){

    ## This function is very similar to 'sim_BM_budding_exp'. The only difference is that we are not repeating the simulation of the cladogenetic events.

    ## Define the size of the chunks to be used for the simulation.
    ## Max age helps to deal with non-ultrametric trees.
    chunk_length <- max( diag( vcv.phylo(sim_BM$contsim) ) ) / 1000 ## 1000 pieces of tree age.

    ## Check if the jump_size is valid.
    if( is.numeric(jump_size) & jump_size < 0.0 ) stop("jump_size cannot be a negative number.")
    if( !is.logical(jump_all) ) stop("jump_all needs to be TRUE or FALSE.")
    jump_type <- match.arg(jump_type, choices=c("fixed", "gaussian"), several.ok=FALSE)

    tt <- reorder.phylo(sim_BM$contsim) ## reorder the tree cladewise.
    ## Need a 'maps' list in the same format as 'phytools' 'simmap' object.
    maps <- vector(mode = "list", length = nrow(tt$edge))

    ## For each edge of the tree we need an starting and an ending state.
    STATES <- matrix(NA, nrow(tt$edge), 2)
    ## For each edge of the tree we need a number to keep track of the lineage.
    ANCESTRY <- vector(mode = "numeric", length = nrow(tt$edge) )

    ## We need to track history instead of creating a new one. So here we define DESTINY:
    DESTINY <- sim_BM$ancestry ## The object from the previous simulation.

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
    STATES[root_id, 1] <- as.numeric( anc )

    ## Initiate the ancestry of the first lineage.
    SCALER[root_id] <- 1.0
    start_scaler <- 1.0

    ## Create a matrix for the storage of the rate scalers associated with the age.
    mat_scaler <- matrix(data = NA, nrow = 1, ncol = 4)

    for (j in 1:nrow(tt$edge)){

        ## ###################################
        ## SIMULATION ALONG THE BRANCH - ANAGENETIC CHANGES.
        ## This need to be done by the chunk.
        chunk_vec <- get_chunk_vec(ll = tt$edge.length[j], chunk_length = chunk_length)
        chunk_state <- vector(mode = "character", length = length(chunk_vec))
        chunk_scaler <- vector(mode = "numeric", length = length(chunk_vec))
        start_state <- as.numeric(STATES[j,1]) ## Continuous value.

        if( sum(ANCESTRY == ANCESTRY[j]) == 1 ){
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
            ## Scale the BM rate and simulate the trait evolution at the chunk.
            sigma_tmp <- sigma * chunk_vec[w] * age_tmp_scaler
            new_state <- as.numeric(start_state) + rnorm(n = 1, sd = sqrt(sigma_tmp) )
            ## Log the state and the scaler values.
            mat_scaler <- rbind(mat_scaler, c(new_lineage, age_tmp, chunk_vec[w], age_tmp_scaler))
            chunk_state[w] <- new_state
            start_state <- new_state ## Update the state at the chunk.
            chunk_scaler[w] <- age_tmp_scaler
        }

        ## Update the maps for the state and the scaler values.
        ## Note that the state here is a numeric value.
        SCALER[[j]] <- data.frame(chunk = chunk_vec, scaler = chunk_scaler)
        maps[[j]] <- setNames(object = chunk_vec, nm = chunk_state)
        new <- as.numeric( chunk_state[ length(chunk_state) ] ) ## The state at the end of the edge.
        ## ###################################

        ## ###################################
        ## STEPS AT THE DESCENDANT NODE - BUDDING SPECIATION.
        ## The state that was sampled now is the state at the descendant node.
        STATES[j, 2] <- as.numeric( new )
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
                ## This is the most important trait of a mother lineage.
                STATES[mother_lineage_id, 1] <- new

                ## Check if we should make a jump at the new species.
                if( is.na(jump_size) ){
                    STATES[daughter_lineage_id, 1] <- new
                } else{
                    ## Make the cladogenetic jump.
                    if( jump_type == "fixed" ){
                        STATES[daughter_lineage_id, 1] <- new + (sample(x = c(1.0,-1.0), size = 1) * new * jump_size)
                    } else if( jump_type == "gaussian" ){
                        STATES[daughter_lineage_id, 1] <- rnorm(n = 1, mean = new, sd = jump_size)
                    } else{
                        stop("Jump type needs to be either fixed or Gaussian.")
                    }
                }

            } else{
                ## Not a budding event, two new lineages will be produced.
                if( jump_all ){
                    ## All new lineages have a jump.
                    if( jump_type == "fixed" ){
                        direction_fix_jump <- sample(x = c(1.0,-1.0), size = length(ii), replace = TRUE)
                        STATES[ii, 1] <- new + (direction_fix_jump * new * jump_size)
                    } else if( jump_type == "gaussian" ){
                        STATES[ii, 1] <- rnorm(n = length(ii), mean = new, sd = jump_size)
                    } else{
                        stop("Jump type needs to be either fixed or Gaussian.")
                    }
                    ANCESTRY[ii] <- DESTINY[ii]
                } else{
                ## Both descendants are new lineages and they share the same state.
                STATES[ii, 1] <- new
                ANCESTRY[ii] <- DESTINY[ii]
                }
            }

        }
        ## ###################################
    }

    ## Get the states for the tips:
    xx_temp <- sapply(1:Ntip(tt), function(n, S, E) S[which(E == n)], S = STATES[, 2], E = tt$edge[,2])
    x <- setNames(object = as.numeric(xx_temp), nm = tt$tip.label)

    ## Temporary object to return the rate scalers.
    colnames( mat_scaler ) <- c("new_lineage", "age", "chunk_length", "scaler")

    return( list(contsim = tt, tip_state = x, edge_state = STATES
                 , ancestry = ANCESTRY, rate_scaler_maps = SCALER
                 , scaler_mat = mat_scaler, maps = maps) )
}

