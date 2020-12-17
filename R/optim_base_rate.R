#' Adjust average rate under budding speciation
#'
#' Budding speciation combined with lineage-age dependent trait evolution will cause the rates of trait evolution to vary across the branches of the tree. This function allows one to adjust the base rate used to simulate traits given a target average rate. This function requires the previous usage of ‘sim_Mk_budding_exp’ or ‘sim_Mk_trace_history’. It is useful to create multiple scenarios of trait evolution while making sure that the average rate across the branch of the tree remains the same.
#' @param simMK the output from ‘sim_Mk_budding_exp’ or ‘sim_Mk_trace_history’
#' @param target_mean_rate a numeric value. The target average rate of trait evolution. At the moment, this function only allows for equal rates (ER) models.
#'
#' @return The base rate parameter which will produce the target average rate.
#' @export
#'
#' @importFrom stats optim
optim_base_rate <- function(simMK, target_mean_rate){
    ## Initial state for the search
    init_par <- runif(n = 1, min = 0.00001, max = 10)
    ## Extract the history from the simulation
    hist_mat <- simMK$scaler_mat[-1,]
    ## The function to be optimized
    dist_fn <- function(par){
        wgt_avg <- apply(X = hist_mat, MARGIN = 1, FUN = function(x) x["chunk_length"] * x["scaler"] * par)
        obs_mean_rate <- sum( wgt_avg ) / sum(hist_mat[,"chunk_length"])
        dist <- sqrt( (target_mean_rate - obs_mean_rate)^2 )
        return( dist )
    }
    ## Run a search
    res <- optim(par = init_par, fn = dist_fn, gr = NULL, method = "Brent", lower = 0.0000001, upper = 100.0)
    return( res )
}

