#' Plot trajectory of continuous trait simulation
#'
#' Function plots the trajectory of a single continuous trait evolving along the branches of the tree following a lineage-age dependent process including events of budding (or symmetric) speciation.
#' @param sim_obj the output of a 'buddPhy' simulation function.
#' @param edge_width the width of the lines
#' @param bg the color for the background of the plot. Helps to visualize some colors. Use "white" to turn it off.
#' @param jump_col the color of the line connecting the states for the nodes in which a phenotype jump happened after budding speciation.
#' @param jump_width the line width for the jump lines.
#' @param main a title for the plot.
#' @param legend if to add a legend. Note that legend can get in the way of the plot sometimes.
#'
#' @return
#' @export
#' @importFrom viridisLite cividis
plot_cont_sim <- function(sim_obj, edge_width = 1.5, bg = "white", jump_col = "coral", jump_width = 1.0, main = "", legend = TRUE){
    ## Simple function to plot the trajectory of the continuous simulation.
    plot_map <- list()
    for( i in 1:length(sim_obj$maps) ){
        ## Can use cumsum to add the time chunks.
        plot_map[[i]] <- c(setNames(object = 0.0, nm = sim_obj$edge_state[i,1]), cumsum(sim_obj$maps[[i]]))
    }
    node_age_order <- node.depth.edgelength(phy = sim_obj$contsim)
    node_age <- vector(mode = "numeric", length = nrow(sim_obj$edge_state) )
    for( i in 1:length(node_age) ){
        nd <- sim_obj$contsim$edge[i,1]
        node_age[i] <- node_age_order[nd]
    }
    range_traits <- range(as.numeric(names(do.call(c, plot_map))))
    range_age <- range( node_age )

    plot(x = 1, xlim = range_age, ylim = range_traits, type = "n", xlab = "Time"
         , ylab = "Trait", main = main)
    rect(xleft = range_age[1], xright = range_age[2], ybottom = range_traits[1]
         , ytop = range_traits[2], col = bg, border = bg)
    if( legend ){
        ## Optional because it can get in the way.
        legend(x = 0.0, y = max( range_traits )
               , legend = seq(from = 0, to = 1, length.out = 11)[-1]
               , fill = cividis(n = 10, direction = -1)
               , ncol = 3)
    }

    ## Prepare the vector in case of jump of character states.
    ## This is the first part to be plotted, so the jump lines stay in the background.
    edge <- sim_obj$contsim$edge
    edge_st <- sim_obj$edge_state
    node_vec <- seq(from = Ntip(sim_obj$contsim)+1, by = 1, length.out = Nnode(sim_obj$contsim))
    jump_id <- sapply(node_vec, function(x) diff(edge_st[which( edge[,1] == x ),1],) != 0 )
    if( any( jump_id ) ){ ## If a jump happened, then plot the lines.
        jump_node <- node_vec[ jump_id ]
        jump_mat <- sapply(jump_node, function(x) edge_st[which( edge[,1] ==  x),1])
        internal_node_age <- node_age_order[ node_vec ] ## Node age also include tips!
        jump_age <- internal_node_age[ jump_id ]
        for( i in 1:length(jump_node) ){
            lines(x = rep(jump_age[i], 2), y = jump_mat[,i], col = jump_col, lwd = jump_width)
        }
    }

    ## Make the plot for the lineages.
    ht_col <- cividis(n = 10, direction = -1)
    for( i in 1:length(plot_map) ){
        trait <-  as.numeric(names(plot_map[[i]]))
        age <- unname(plot_map[[i]]) + node_age[i]
        line_levels <- cut(x = sim_obj$rate_scaler_maps[[i]]$scaler
                           , breaks = seq(from = 0, to = 1, length.out = 11)
                           , labels = FALSE)
        line_cols <- ht_col[line_levels]
        for( i in seq(from = 1, to = length(trait)-1) ){
            lines(y = trait[i:(i+1)], x = age[i:(i+1)], col = line_cols[i], lwd = edge_width)
        }
    }

}
