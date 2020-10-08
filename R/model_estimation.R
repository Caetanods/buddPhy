#' Rescale branch lengths following the model
#'
#' Function to rescale branches of a phylogenetic tree following the exponential rate proportional, or inversely-proportional, to lineage-age. This function can be used to simulate a realization of the model for posterior predictive checks among other reasons.
#' @param change_rate rate for the exponential distribution.
#' @param budding_prob probability of a budding speciation on the internal node.
#' @param tree phylogenetic tree.
#' @param chunk_fraction fraction of the total tree-height for branch length discretization.
#' @param decay_fn whether rates should vary proportional (FALSE) or inversely proportional (TRUE) to lineage age.
#'
#' @return a phylogenetic tree with rescaled branch lengths and reordered to "cladewise" order (see 'ape::reorder.phylo').
#' @export
tree_rescale_exp <- function(change_rate = 0.05, budding_prob = 0.3, tree, chunk_fraction = 0.01, decay_fn = TRUE){
    ## FREE PARAMETERS: change_rate and budding_prob.
    ## tree, chunk_fraction, and decay_fn are constants in the model.
    ## Produces a data augmentation map for both the rate scalers and the mother lineages across the tree.

    ## Define the size of the chunks to be used for the simulation.
    chunk_length <- vcv.phylo(tree)[1,1] * chunk_fraction

    tt <- reorder.phylo(x = phy, order = "cladewise") ## reorder the tree cladewise.
    rescaled_tree <- tt ## Tree to be rescaled.
    root_id <- which(tt$edge[, 1] == Ntip(tt) + 1) ## The root edges.

    ## For each edge of the tree we need a number to keep track of the lineage.
    ANCESTRY <- vector(mode = "numeric", length = nrow(tt$edge) )
    ## Initiate the ancestry of the first lineage. At the edges connecting to the root we have the first and the second lineages
    ANCESTRY[root_id] <- c(1,2)
    scaler_edge_vec <- vector(mode = "numeric", length = nrow(tt$edge))
    scaler_edge_vec[root_id] <- 1.0 ## Starting scaler for either decay or positive lineage-age rates models.
    ## This is an homogeneous decay process, the starting scaler is always the same.
    start_scaler <- 1.0

    for (j in 1:nrow(tt$edge)){

        ## For each edge we use the scaling function to rescale the branch.
        chunk_vec <- get_chunk_vec(ll = tt$edge.length[j], chunk_length = chunk_length)
        rescaled_chunk_vec <- chunk_vec ## Chunks to be rescaled.

        if( sum(ANCESTRY == ANCESTRY[j]) == 1 ){
            ## First time the lineage number appears in the ANCESTRY vector.
            ## New lineage, including the root edges.
            ## start_scaler <- 1.0
            new_lineage <- TRUE
        } else{
            ## Mother lineage. Lineage already exists.
            ## The ancestral edge NEEDS to belong to the same continued lineage.
            ancestral_edge <- tt$edge[,2] == tt$edge[j,1]
            ## Same as line below, but here we just want the value at the end of the branch.
            ## start_scaler <- scaler_edge_vec[ancestral_edge]
            ## start_scaler <- get_des_scaler(scaler = SCALER, id = ancestral_edge)
            new_lineage <- FALSE
        }

        for( w in 1:length(chunk_vec) ){
            ## Implement the lineage age function:
            chunk_sum <- ifelse(test = w == 1, yes = 0, no = sum(chunk_vec[1:(w-1)]) )
            ## get_chunk_age will compute the age of the current chunk.
            age_tmp <- get_chunk_age(tree = tt, current_node = tt$edge[j,1], chunk_length = chunk_vec[w], chunk_sum = chunk_sum, lineage_number = ANCESTRY[j], ancestry = ANCESTRY)
            ## Take one step on the scaler based on the change_rate and the value of the previous scaler.
            ## If the change_rate is 0.0, then this scaler will always be 1.0 and the model will be time-homogeneous.
            if( decay_fn ){
                ## The relative start scaler never changes.
                chunk_scaler <- start_scaler / exp(change_rate * age_tmp)
            } else{
                chunk_scaler <- start_scaler * exp(change_rate * age_tmp)
            }
            ## Scale the chunk of the branch:
            ## This will need to be stored in a copy of the tree, otherwise we will get wrong chunk ages.
            rescaled_chunk_vec[w] <- rescaled_chunk_vec[w] * chunk_scaler
        }

        ## Log the scaler value at the end of edge j.
        ## This substitutes the SCALER list.
        scaler_edge_vec[j] <- chunk_scaler

        ## Update the branch of the rescaled tree. This will be the sum of the chunks.
        rescaled_tree$edge.length[j] <- sum( rescaled_chunk_vec )

        ## ###################################
        ## STEPS AT THE DESCENDANT NODE - BUDDING SPECIATION.
        ## Check if the descendant node is an internal node or tip.
        ii <- which(tt$edge[, 1] == tt$edge[j, 2])
        if(length(ii) > 0){ ## Internal node (two descendant edges).
            ## Check if budding speciation happens.
            budd_happens <- sample(x = c(TRUE, FALSE), size = 1, prob = c(budding_prob, 1-budding_prob) )

            ## #####################################
            ## Apply budding speciation if necessary.
            if( budd_happens ){
                ## One lineage keeps the same ancestry (at random).
                ndrd <- sample(x = c(1,2), size = 2, replace = FALSE)
                ## Update the continuation of the ancestral.
                ANCESTRY[ii[ndrd[1]]] <- ANCESTRY[j] ## j is the ancestral edge.
                ## Update the new lineage, with the cladogenetic function.
                ANCESTRY[ii[ndrd[2]]] <- max(ANCESTRY) + 1 ## the new lineage.
            } else{
                ## Both descendants are new lineages and they share the same state.
                ANCESTRY[ii] <- c(1,2) + max( ANCESTRY )
            }
        }
    }

    ## Return the rescaled tree:
    return( rescaled_tree )
}

##' @noRd
initialize_budding_history <- function(tree, change_rate = 0.05, budding_prob = 0.3, chunk_fraction = 0.01, decay_fn = TRUE){
    ## Generates a random budding history from the prior.
    ## Returns budding nodes,ancestry, and the rescaled tree.

    ## FREE PARAMETERS: change_rate and budding_prob.
    ## tree, chunk_fraction, and decay_fn are constants in the model.

    ## Define the size of the chunks to be used for the simulation.
    chunk_length <- vcv.phylo(tree)[1,1] * chunk_fraction

    tt <- reorder.phylo(x = tree, order = "cladewise") ## reorder the tree cladewise.
    rescaled_tree <- tt ## Tree to be rescaled.
    root_id <- which(tt$edge[, 1] == Ntip(tt) + 1) ## The root edges.

    ## For each edge of the tree we need a number to keep track of the lineage.
    ANCESTRY <- vector(mode = "numeric", length = nrow(tt$edge) )
    ## For each edge of the tree we keep if the node is or not a budding speciation type.
    BUDD <- rep(x = FALSE, times = nrow(tt$edge) )
    ## Initiate the ancestry of the first lineage. At the edges connecting to the root we have the first and the second lineages
    ANCESTRY[root_id] <- c(1,2)
    scaler_edge_vec <- vector(mode = "numeric", length = nrow(tt$edge))
    scaler_edge_vec[root_id] <- 1.0 ## Starting scaler for either decay or positive lineage-age rates models.
    ## Starting scaler is always the same. Homogeneos process.
    start_scaler <- 1.0

    for (j in 1:nrow(tt$edge)){
        ## For each edge we use the scaling function to rescale the branch.
        chunk_vec <- buddPhy:::get_chunk_vec(ll = tt$edge.length[j], chunk_length = chunk_length)
        rescaled_chunk_vec <- chunk_vec ## Chunks to be rescaled.

        if( sum(ANCESTRY == ANCESTRY[j]) == 1 ){
            ## First time the lineage number appears in the ANCESTRY vector.
            ## New lineage, including the root edges.
            ## start_scaler <- 1.0
            new_lineage <- TRUE
        } else{
            ## Mother lineage. Lineage already exists.
            ## The ancestral edge NEEDS to belong to the same continued lineage.
            ancestral_edge <- tt$edge[,2] == tt$edge[j,1]
            ## Same as line below, but here we just want the value at the end of the branch.
            ## start_scaler <- scaler_edge_vec[ancestral_edge]
            ## start_scaler <- get_des_scaler(scaler = SCALER, id = ancestral_edge)
            new_lineage <- FALSE
        }

        for( w in 1:length(chunk_vec) ){
            ## Implement the lineage age function:
            chunk_sum <- ifelse(test = w == 1, yes = 0, no = sum(chunk_vec[1:(w-1)]) )
            ## get_chunk_age will compute the age of the current chunk.
            age_tmp <- buddPhy:::get_chunk_age(tree = tt, current_node = tt$edge[j,1], chunk_length = chunk_vec[w], chunk_sum = chunk_sum, lineage_number = ANCESTRY[j], ancestry = ANCESTRY)
            ## Take one step on the scaler based on the change_rate and the value of the previous scaler.
            ## If the change_rate is 0.0, then this scaler will always be 1.0 and the model will be time-homogeneous.
            if( decay_fn ){
                chunk_scaler <- start_scaler / exp(change_rate * age_tmp)
            } else{
                chunk_scaler <- start_scaler * exp(change_rate * age_tmp)
            }
            ## Scale the chunk of the branch:
            ## This will need to be stored in a copy of the tree, otherwise we will get wrong chunk ages.
            rescaled_chunk_vec[w] <- rescaled_chunk_vec[w] * chunk_scaler
        }

        ## Log the scaler value at the end of edge j.
        ## This substitutes the SCALER list.
        scaler_edge_vec[j] <- chunk_scaler

        ## Update the branch of the rescaled tree. This will be the sum of the chunks.
        rescaled_tree$edge.length[j] <- sum( rescaled_chunk_vec )

        ## STEPS AT THE DESCENDANT NODE - BUDDING SPECIATION.
        ## Check if the descendant node is an internal node or tip.
        ii <- which(tt$edge[, 1] == tt$edge[j, 2])
        if(length(ii) > 0){ ## Internal node (two descendant edges).
            ## Check if budding speciation happens.
            BUDD[j] <- sample(x = c(TRUE, FALSE), size = 1, prob = c(budding_prob, 1-budding_prob) )

            ## #####################################
            ## Apply budding speciation if necessary.
            if( BUDD[j] ){
                ## One lineage keeps the same ancestry (at random).
                ndrd <- sample(x = c(1,2), size = 2, replace = FALSE)
                ## Update the continuation of the ancestral.
                ANCESTRY[ii[ndrd[1]]] <- ANCESTRY[j] ## j is the ancestral edge.
                ## Update the new lineage, with the cladogenetic function.
                ANCESTRY[ii[ndrd[2]]] <- max(ANCESTRY) + 1 ## the new lineage.
            } else{
                ## Both descendants are new lineages and they share the same state.
                ANCESTRY[ii] <- c(1,2) + max( ANCESTRY )
            }
        }
    }

    ## Return the objects needed.
    return( list(res_tree = rescaled_tree, ANCESTRY = ANCESTRY, BUDD = BUDD, SCALER = scaler_edge_vec ) )
}

##' @noRd
update_budding_history <- function(budd_hist, change_rate = 0.05, budding_prob = 0.3, tree, chunk_fraction = 0.01, decay_fn = TRUE, skip_shallow = TRUE){
    ## Choose a node, updates the history of the lineages descending this node.
    ## The other strategy would be to modify the current history. However, because all lineages are connected and the mother lineage history depends on the continuation of these lineages, it will be complicated, and arbitrary, the way that we modify them.
    ## Simulating the history again, following the budding_prob parameter, seems a better strategy.
    ## budd_hist is a list following the output from initialize_budding_history .
    ## Need all the parameters to be able to return another rescaled tree following the updated history.

    ## Recover the current lineage history.
    ANCESTRY <- budd_hist$ANCESTRY
    BUDD <- budd_hist$BUDD
    scaler_edge_vec <- budd_hist$SCALER
    edge <- budd_hist$res_tree$edge
    tree <- budd_hist$res_tree

    ## Sample node to change:
    nodes_vec <- seq(from = Ntip(tree)+1, length.out = Nnode(tree))
    if( skip_shallow ){
        ## Skip the shallow nodes which BOTH descendants are tips.
        drop_nodes <- sapply(nodes_vec, function(x) all( edge[ edge[,1] == x, 2] <= Ntip(tree) ) )
        nodes_vec <- nodes_vec[ !drop_nodes ]
    }

    prop_nd <- sample(x = nodes_vec, size = 1)

    ## Find all the edges descending the sampled node.
    des_nodes <- phytools::getDescendants(tree = tree, node = prop_nd)
    des_nodes <- des_nodes[ des_nodes > Ntip(tree) ]
    row_update <- which( edge[,1] %in% c(des_nodes, prop_nd) )

    ## Erase the history AFTER the first chosen node.
    ## ANCESTRY tracks the stem of the nodes. Each node is represented twice in this vector.
    reset_anc <- row_update[ edge[row_update,1] != prop_nd ]
    ANCESTRY[ reset_anc ] <- 0

    ## ########################################################
    ## Getting some extra objects for the simulation.
    chunk_length <- vcv.phylo(tree)[1,1] * chunk_fraction
    tt <- reorder.phylo(tree) ## reorder the tree cladewise.
    rescaled_tree <- tt ## Tree to be rescaled.
    start_scaler <- 1.0 ## Start scaler is always the same.

    ## Update all descending lineages:
    ## for (j in 1:nrow(tt$edge)){
    for(j in row_update){

        ## For each edge we use the scaling function to rescale the branch.
        chunk_vec <- buddPhy:::get_chunk_vec(ll = tt$edge.length[j], chunk_length = chunk_length)
        rescaled_chunk_vec <- chunk_vec ## Chunks to be rescaled.

        if( sum(ANCESTRY == ANCESTRY[j]) == 1 ){
            ## First time the lineage number appears in the ANCESTRY vector.
            ## New lineage, including the root edges.
            ## start_scaler <- 1.0
            new_lineage <- TRUE
        } else{
            ## Mother lineage. Lineage already exists.
            ## The ancestral edge NEEDS to belong to the same continued lineage.
            ancestral_edge <- tt$edge[,2] == tt$edge[j,1]
            ## Same as line below, but here we just want the value at the end of the branch.
            ## start_scaler <- scaler_edge_vec[ancestral_edge]
            ## start_scaler <- get_des_scaler(scaler = SCALER, id = ancestral_edge)
            new_lineage <- FALSE
        }

        for( w in 1:length(chunk_vec) ){
            ## Implement the lineage age function:
            chunk_sum <- ifelse(test = w == 1, yes = 0, no = sum(chunk_vec[1:(w-1)]) )
            ## get_chunk_age will compute the age of the current chunk.
            age_tmp <- buddPhy:::get_chunk_age(tree = tt, current_node = tt$edge[j,1], chunk_length = chunk_vec[w], chunk_sum = chunk_sum, lineage_number = ANCESTRY[j], ancestry = ANCESTRY)
            ## Take one step on the scaler based on the change_rate and the value of the previous scaler.
            ## If the change_rate is 0.0, then this scaler will always be 1.0 and the model will be time-homogeneous.
            if( decay_fn ){
                chunk_scaler <- start_scaler / exp(change_rate * age_tmp)
            } else{
                chunk_scaler <- start_scaler * exp(change_rate * age_tmp)
            }
            ## Scale the chunk of the branch:
            ## This will need to be stored in a copy of the tree, otherwise we will get wrong chunk ages.
            rescaled_chunk_vec[w] <- rescaled_chunk_vec[w] * chunk_scaler
        }

        ## Log the scaler value at the end of edge j.
        ## This substitutes the SCALER list.
        scaler_edge_vec[j] <- chunk_scaler

        ## Update the branch of the rescaled tree. This will be the sum of the chunks.
        rescaled_tree$edge.length[j] <- sum( rescaled_chunk_vec )

        ## STEPS AT THE DESCENDANT NODE - BUDDING SPECIATION.
        ## Check if the descendant node is an internal node or tip.
        ii <- which(tt$edge[, 1] == tt$edge[j, 2])
        if(length(ii) > 0){ ## Internal node (two descendant edges).
            ## Check if budding speciation happens.
            BUDD[j] <- sample(x = c(TRUE, FALSE), size = 1, prob = c(budding_prob, 1-budding_prob) )

            ## #####################################
            ## Apply budding speciation if necessary.
            if( BUDD[j] ){
                ## One lineage keeps the same ancestry (at random).
                ndrd <- sample(x = c(1,2), size = 2, replace = FALSE)
                ## Update the continuation of the ancestral.
                ANCESTRY[ii[ndrd[1]]] <- ANCESTRY[j] ## j is the ancestral edge.
                ## Update the new lineage, with the cladogenetic function.
                ANCESTRY[ii[ndrd[2]]] <- max(ANCESTRY) + 1 ## the new lineage.
            } else{
                ## Both descendants are new lineages and they share the same state.
                ANCESTRY[ii] <- c(1,2) + max( ANCESTRY )
            }
        }
    }

    ## Return the objects needed.
    return( list(res_tree = rescaled_tree, ANCESTRY = ANCESTRY, BUDD = BUDD, SCALER = scaler_edge_vec ) )
}

##' @noRd
update_branch_scaling <- function(budd_hist, tree, change_rate = 0.05, chunk_fraction = 0.01, decay_fn = TRUE){
    ## budd_hist is a list following the output from initialize_budding_history
    ## Keep the lineage history the same, but rescale the branch lengths.

    ## Define the size of the chunks to be used for the simulation.
    chunk_length <- vcv.phylo(tree)[1,1] * chunk_fraction

    tt <- reorder.phylo(tree) ## reorder the tree cladewise.
    rescaled_tree <- tt ## Tree to be rescaled.
    root_id <- which(tt$edge[, 1] == Ntip(tt) + 1) ## The root edges.

    ## The lineage history will keep constant.
    ## These objects have the "future" of the simulation, but we need to ge step by step.
    ## [[ Otherwise I will need to change the code too much. This is a hack. ]]
    real_ANCESTRY <- budd_hist$ANCESTRY
    real_BUDD <- budd_hist$BUDD
    ## The "observed" vectors for the for loop.
    ANCESTRY <- vector(mode = "numeric", length = nrow(tt$edge) )
    ANCESTRY[root_id] <- c(1,2)
    BUDD <- rep(x = FALSE, times = nrow(tt$edge) )

    ## We will sample the scalers.
    scaler_edge_vec <- vector(mode = "numeric", length = nrow(tt$edge))
    scaler_edge_vec[root_id] <- 1.0 ## Starting scaler for either decay or positive lineage-age rates models.
    start_scaler <- 1.0

    for(j in 1:nrow(tt$edge)){

        ## For each edge we use the scaling function to rescale the branch.
        chunk_vec <- buddPhy:::get_chunk_vec(ll = tt$edge.length[j], chunk_length = chunk_length)
        rescaled_chunk_vec <- chunk_vec ## Chunks to be rescaled.

        if( sum(ANCESTRY == ANCESTRY[j]) == 1 ){
            ## First time the lineage number appears in the ANCESTRY vector.
            ## New lineage, including the root edges.
            ## start_scaler <- 1.0
            new_lineage <- TRUE
        } else{
            ## Mother lineage. Lineage already exists.
            ## The ancestral edge NEEDS to belong to the same continued lineage.
            ancestral_edge <- tt$edge[,2] == tt$edge[j,1]
            ## Same as line below, but here we just want the value at the end of the branch.
            ## start_scaler <- scaler_edge_vec[ancestral_edge]
            ## start_scaler <- get_des_scaler(scaler = SCALER, id = ancestral_edge)
            new_lineage <- FALSE
        }

        for( w in 1:length(chunk_vec) ){
            ## Implement the lineage age function:
            chunk_sum <- ifelse(test = w == 1, yes = 0, no = sum(chunk_vec[1:(w-1)]) )
            ## get_chunk_age will compute the age of the current chunk.
            age_tmp <- buddPhy:::get_chunk_age(tree = tt, current_node = tt$edge[j,1], chunk_length = chunk_vec[w], chunk_sum = chunk_sum, lineage_number = ANCESTRY[j], ancestry = ANCESTRY)
            ## Take one step on the scaler based on the change_rate and the value of the previous scaler.
            ## If the change_rate is 0.0, then this scaler will always be 1.0 and the model will be time-homogeneous.
            if( decay_fn ){
                chunk_scaler <- start_scaler / exp(change_rate * age_tmp)
            } else{
                chunk_scaler <- start_scaler * exp(change_rate * age_tmp)
            }
            ## Scale the chunk of the branch:
            ## This will need to be stored in a copy of the tree, otherwise we will get wrong chunk ages.
            rescaled_chunk_vec[w] <- rescaled_chunk_vec[w] * chunk_scaler
        }

        ## Log the scaler value at the end of edge j.
        ## This substitutes the SCALER list.
        scaler_edge_vec[j] <- chunk_scaler

        ## Update the branch of the rescaled tree. This will be the sum of the chunks.
        rescaled_tree$edge.length[j] <- sum( rescaled_chunk_vec )

        ## Update the observed BUDD and ANCESTRY objects.
        ii <- which(tt$edge[, 1] == tt$edge[j, 2])
        BUDD[j] <- real_BUDD[j]
        ANCESTRY[ii] <- real_ANCESTRY[ii]
    }

    ## Return the objects needed.
    return( list(res_tree = rescaled_tree, ANCESTRY = ANCESTRY, BUDD = BUDD, SCALER = scaler_edge_vec ) )
}

##' @noRd
re_build_Q <- function( q_vec, n_states, phytools = FALSE){
    ## Helping function to transform the Q matrix back and forth.
    build_Q <- matrix(nrow = n_states, ncol = n_states)
    ii <- 1
    for( i in 1:n_states ){
        for( j in 1:n_states ){
            if( i != j ){
                if( phytools ){
                    build_Q[j,i] <- q_vec[ii]
                } else{
                    build_Q[i,j] <- q_vec[ii]
                }
                ii <- ii + 1
            }
        }
    }
    diag( build_Q ) <- -1 * rowSums( build_Q, na.rm = TRUE)
    return( build_Q )
}

##' @noRd
get_edge <- function(hh){
    ## Need to transform the tree before getting the edge lengths to compute the likelihood.
    ## This function will help with the transformation.
    phy <- hh$res_tree
    phy_ord <- reorder.phylo(x = phy, order = "postorder")
    edge <- phy_ord$edge.length
    return( edge )
}

#' Bayesian estimation of discrete trait evolution with budding speciation
#'
#' Estimates rates of discrete trait evolution under a Markov model with a single transition matrix fitted to the tree. Implements a lineage-age model which influences how fast trait evolution happens given the age of each lineage. Lineage-age is controlled by a budding probability factor. Budding speciation, as applied here, means that one of the descendant lineages from an internal node will be a continuation of the ancestral lineage. The other descendant will be a new lineage.
#'
#' @param phy phylogeny with branch lengths proportional to time
#' @param trait named vector with discrete traits
#' @param k a character vector with the names of the states. Need to match the states in 'traits'
#' @param prior_beta_Q a numeric vector with two elements informing the alpha and beta shape parameters for a Beta distribution used as a prior to the rates of transition. Here we scale the Beta distribution such that the maximum rate will be equal to 10x the transition rate estimated for a all equal rates model, and used as the starting point for the MCMC estimation.
#' @param prior_beta_budd a numeric vector with two elements informing the alpha and beta shape parameters for a Beta distribution used as a prior to the probability of budding speciation applied at each node of the phylogeny.
#' @param prior_exp_rate the mean of the exponential distribution used as a prior for the exponential decay rate applied to the model.
#' @param sample_prob a numeric vector with the proposal ratio for each of the four free parameters of the model. In order: transition rates, budding probabilities, exponential rate function, lineage history.
#' @param Q_prop proposal rate for a multiplier rate proposal
#' @param budd_prop proposal rate for a multiplier rate proposal
#' @param rate_prop proposal rate for a multiplier rate proposal
#' @param gen total number of generations, including burn-in, used in the MCMC estimation
#'
#' @return a list with the acceptance rate (accept), the likelihood (lik), the prior probabilities (prior), the transition rates (Q), the probability of budding speciation(budd), and the posterior distribution of the rate scalers (exp_rate).
#' @export
#'
mcmc_mk_budd_exp <- function(phy, trait, k, prior_beta_Q = c(2,5), prior_beta_budd = c(2,5), prior_exp_rate = 1, sample_prob = c(0.25,0.25,0.25,0.25), Q_prop = 0.2, budd_prop = 0.2, rate_prop = 0.2, gen = 100){
    ## Order of sampling: c("Q","budd_prob","rate_fn","history")

    ## Make preliminary checks:
    sample_prob <- sample_prob / sum(sample_prob)

    ## Prepare constants to run ratematrix likelihood function.
    phy <- ape::reorder.phylo(x = phy, order = "postorder")
    n_nodes <- Nnode(phy)
    n_states <- length(k)
    n_tips <- Ntip(phy)
    edge_mat <- phy$edge
    parents <- unique( phy$edge[,1] )
    ## "X" makes sure that the Q matrix is in the same order of the states "k".
    X <- ratematrix:::makeDataTips(X = trait, states = k)
    root_node <- n_tips + 1
    root_type <- 1

    ## Define containers ####
    budd_chain <- numeric(length = gen) ## Budding probability.
    exp_rate_chain <- numeric(length = gen) ## Exponential parameter.
    k_matrix <- (n_states * n_states) - n_states ## The elements on the diversitree likelihood.
    Q_chain <- matrix(data = NA, nrow = gen, ncol = k_matrix) ## Transition matrix parameters.
    edge_chain <- matrix(data = NA, nrow = gen, ncol = length(phy$edge.length))
    log_lik_chain <- numeric(length = gen) ## log likelihood.
    log_prior_chain <- numeric(length = gen) ## log prior.
    accept_chain <- numeric(length = gen) ## Acceptance.
    ## Add a container for the budd nodes and lineage history as well.

    ## Starting point ####
    ## Using MLE to get a quick estimate of the rates.
    start_mk <- phytools::fitMk(tree = phy, x = trait, model = "ER")
    Q_chain[1,] <- rep(start_mk$rates, times = k_matrix)
    budd_chain[1] <- rbeta(n = 1, shape1 = prior_beta_budd[1], shape2 = prior_beta_budd[2])
    exp_rate_chain[1] <- rexp(n = 1, rate = prior_exp_rate)
    accept_chain[1] <- 1

    ## Get the starting point for the edge.length. Need to transform the tree.
    ## Here we also get the current state of the lineage history.
    hist_curr <- buddPhy:::initialize_budding_history(tree = phy, change_rate = exp_rate_chain[1]
                                                      , budding_prob = budd_chain[1]
                                                      , chunk_fraction = 0.01, decay_fn = TRUE)
    edge_chain[1,] <- buddPhy:::get_edge(hh = hist_curr)

    ## The prior probability for the Q rates is a scaled beta distribution.
    q_scale <- start_mk$rates * 10 ## Max value for the rates.
    prior_fn <-function(q_vec, budd, rate){
        ## Function conditioned on constants defined on the environment.
        q_pp <-dbeta(x = q_vec/q_scale, shape1 = prior_beta_Q[1], shape2 = prior_beta_Q[2], log = TRUE)
        budd_pp <- dbeta(x = budd, shape1 = prior_beta_budd[1], shape2 = prior_beta_budd[2], log = TRUE)
        rate_pp <- dexp(x = budd, rate = prior_exp_rate, log = TRUE)
        return( sum(q_pp, budd_pp, rate_pp) )
    }

    ## Init lik and prior.
    ## edge_len and Q are estimated. Need to rebuild Q to get the likelihood.
    log_lik_chain[1] <- ratematrix:::logLikMk_C(n_nodes = n_nodes, n_tips = n_tips
                                                , n_states = n_states
                                                , edge_len = buddPhy:::get_edge(hh = hist_curr)
                                                , edge_mat = edge_mat, parents = parents, X = X
                                                , Q = buddPhy:::re_build_Q(Q_chain[1,], n_states)
                                                , root_node = root_node, root_type = root_type)

    log_prior_chain[1] <- prior_fn(q_vec = Q_chain[1,], budd = budd_chain[1], rate = exp_rate_chain[1])

    ## MCMC loop ####
    for( i in 2:gen ){

        ## Choose the parameter.
        par.update <- sample(x = c("Q","budd_prob","rate_fn","history"), size = 1, prob = sample_prob)
        ## Q - the transition matrix.
        ## budd_prob - the probability of budding events.
        ## rate_fn - the parameter for the exponential decay function.
        ## history - the mother lineages history (the data augmentation part).

        ## Update step ####
        ## We can use 'ratematrix' update functions here too.
        if( par.update == "Q"){
            ## Update a varying number of elements on the Q matrix.
            which_update <- sample(x = 1:k_matrix, size = sample(x = 1:k_matrix, size = 1))
            Q_mult_prop <- sapply(which_update, function(qj)
                ratematrix:::multiplierProposal(x = Q_chain[i-1,qj], a = Q_prop)
            )
            Q_chain_prop <- Q_chain[i-1,]
            Q_chain_prop[which_update] <- Q_mult_prop[1,]

            ## Odds ratio: ####
            lik_curr <- log_lik_chain[i-1]
            lik_prop <- ratematrix:::logLikMk_C(n_nodes = n_nodes, n_tips = n_tips
                                                , n_states = n_states
                                                , edge_len = buddPhy:::get_edge(hh = hist_curr)
                                                , edge_mat = edge_mat, parents = parents, X = X
                                                , Q = buddPhy:::re_build_Q(Q_chain_prop, n_states)
                                                , root_node = root_node, root_type = root_type)
            prior_curr <- log_prior_chain[i-1]
            prior_prop <- prior_fn(q_vec = Q_chain_prop, budd = budd_chain[i-1], rate = exp_rate_chain[i-1])
            ll <- lik_prop - lik_curr
            pp <- prior_prop - prior_curr
            ## ratio of log(lik) + ratio of log(prior) + the proposal ratio (for the multiplier!).
            rr <- ll + pp + sum( Q_mult_prop[2,] )
            ## Accept or reject ####
            ## This is a simple step, because the log will be done after.
            accept_par <- ifelse(test = exp(rr) > runif(1), yes = "Q", no = NA)
        } else if( par.update == "budd_prob"){
            budd_mult_prop <- ratematrix:::multiplierProposal(x = budd_chain[i-1], a = budd_prop)
            budd_update <- budd_mult_prop[1]
            ## Update the entire history given the new value of the budding parameter.
            hist_prop <- buddPhy:::initialize_budding_history(tree = phy, change_rate = exp_rate_chain[i-1]
                                                              , budding_prob = budd_update
                                                              , chunk_fraction = 0.01, decay_fn = TRUE)
            ## Odds ratio: ####
            lik_curr <- log_lik_chain[i-1]
            lik_prop <- ratematrix:::logLikMk_C(n_nodes = n_nodes, n_tips = n_tips, n_states = n_states
                                                , edge_len = buddPhy:::get_edge(hh = hist_prop)
                                                , edge_mat = edge_mat, parents = parents, X = X
                                                , Q = buddPhy:::re_build_Q( Q_chain[i-1,], n_states)
                                                , root_node = root_node, root_type = root_type)
            prior_curr <- log_prior_chain[i-1]
            prior_prop <- prior_fn(q_vec = Q_chain[i-1,], budd = budd_update, rate = exp_rate_chain[i-1])
            ll <- lik_prop - lik_curr
            pp <- prior_prop - prior_curr
            ## ratio of log(lik) + ratio of log(prior) + the proposal ratio (for the multiplier!).
            rr <- ll + pp + budd_mult_prop[2]
            ## Accept or reject ####
            ## This is a simple step, because the log will be done after.
            accept_par <- ifelse(test = exp(rr) > runif(1), yes = "budd_prob", no = NA)
        } else if(par.update == "rate_fn"){
            exp_rate_mult_prop <- ratematrix:::multiplierProposal(x = exp_rate_chain[i-1], a = rate_prop)
            exp_rate_update <- exp_rate_mult_prop[1]
            ## Need to update only the branch scale, keeping the lineage history constant.
            hist_prop <- buddPhy:::update_branch_scaling(budd_hist = hist_curr, tree = phy
                                                         , change_rate = exp_rate_update
                                                         , chunk_fraction = 0.01, decay_fn = TRUE)
            ## Odds ratio: ####
            lik_curr <- log_lik_chain[i-1]
            lik_prop <- ratematrix:::logLikMk_C(n_nodes = n_nodes, n_tips = n_tips, n_states = n_states
                                                , edge_len = buddPhy:::get_edge(hh = hist_prop)
                                                , edge_mat = edge_mat, parents = parents, X = X
                                                , Q = buddPhy:::re_build_Q( Q_chain[i-1,], n_states)
                                                , root_node = root_node, root_type = root_type)
            prior_curr <- log_prior_chain[i-1]
            prior_prop <- prior_fn(q_vec = Q_chain[i-1,], budd = budd_chain[i-1], rate = exp_rate_update)
            ll <- lik_prop - lik_curr
            pp <- prior_prop - prior_curr
            ## ratio of log(lik) + ratio of log(prior) + the proposal ratio (for the multiplier!).
            rr <- ll + pp + exp_rate_mult_prop[2]
            ## Accept or reject ####
            ## This is a simple step, because the log will be done after.
            accept_par <- ifelse(test = exp(rr) > runif(1), yes = "rate_fn", no = NA)
        } else if(par.update == "history"){
            ## Here we are just updating the history by changing the map on some nodes only.
            ## Note that all nodes descending the one chosen are also updated. So many nodes can be updated.
            hist_prop <- buddPhy:::update_budding_history(budd_hist = hist_curr, tree = phy
                                                          , change_rate = exp_rate_chain[i-1]
                                                          , budding_prob = budd_chain[i-1]
                                                          , chunk_fraction = 0.01, decay_fn = TRUE)
            ## Odds ratio: ####
            lik_curr <- log_lik_chain[i-1]
            lik_prop <- ratematrix:::logLikMk_C(n_nodes = n_nodes, n_tips = n_tips, n_states = n_states
                                                , edge_len = buddPhy:::get_edge(hh = hist_prop)
                                                , edge_mat = edge_mat, parents = parents, X = X
                                                , Q = buddPhy:::re_build_Q( Q_chain[i-1,], n_states)
                                                , root_node = root_node, root_type = root_type)
            prior_curr <- log_prior_chain[i-1]
            prior_prop <- prior_fn(q_vec = Q_chain[i-1,], budd = budd_chain[i-1], rate = exp_rate_chain[i-1])
            ll <- lik_prop - lik_curr
            pp <- prior_prop - prior_curr ## pp <- 0.0
            ## ratio of log(lik) + ratio of log(prior)
            rr <- ll + pp
            ## Accept or reject ####
            ## This is a simple step, because the log will be done after.
            accept_par <- ifelse(test = exp(rr) > runif(1), yes = "history", no = NA)
        }

        ## Log the MCMC ####
        ## The log type will depend on the type of accept. The reject is always the same.
        if( is.na( accept_par ) ){
            ## Nothing to update.
            accept_chain[i] <- 0
            log_lik_chain[i] <- log_lik_chain[i-1]
            log_prior_chain[i] <- log_prior_chain[i-1]
            Q_chain[i,] <- Q_chain[i-1,]
            budd_chain[i] <- budd_chain[i-1]
            exp_rate_chain[i] <- exp_rate_chain[i-1]
            ## hist_curr <- hist_curr
        } else if( accept_par == "Q" ){
            ## Update Q
            accept_chain[i] <- 1
            log_lik_chain[i] <- lik_prop
            log_prior_chain[i] <- prior_prop
            Q_chain[i,] <- Q_chain_prop
            budd_chain[i] <- budd_chain[i-1]
            exp_rate_chain[i] <- exp_rate_chain[i-1]
            ## hist_curr <- hist_curr
        } else if( accept_par == "budd_prob" ){
            ## Update budd_prob
            accept_chain[i] <- 1
            log_lik_chain[i] <- lik_prop
            log_prior_chain[i] <- prior_prop
            Q_chain[i,] <- Q_chain[i-1,]
            budd_chain[i] <- budd_update
            exp_rate_chain[i] <- exp_rate_chain[i-1]
            hist_curr <- hist_prop
        } else if( accept_par == "rate_fn" ){
            ## Update rate_fn
            accept_chain[i] <- 1
            log_lik_chain[i] <- lik_prop
            log_prior_chain[i] <- prior_prop
            Q_chain[i,] <- Q_chain[i-1,]
            budd_chain[i] <- budd_chain[i-1]
            exp_rate_chain[i] <- exp_rate_update
            hist_curr <- hist_prop
        } else {
            ## Update history
            accept_chain[i] <- 1
            log_lik_chain[i] <- lik_prop
            log_prior_chain[i] <- prior_prop
            Q_chain[i,] <- Q_chain[i-1,]
            budd_chain[i] <- budd_chain[i-1]
            exp_rate_chain[i] <- exp_rate_chain[i-1]
            hist_curr <- hist_prop
        }

    } ## End of the MCMC loop.

    return( list(accept = accept_chain, lik = log_lik_chain, prior = log_prior_chain, Q = Q_chain
                 , budd = budd_chain, exp_rate = exp_rate_chain) )
}
