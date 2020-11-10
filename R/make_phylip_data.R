#' Creates BioGeoBEARS phylip file
#'
#' @param data a named vector with the tip data.
#' @param states a vector with the states of the data.
#' @param file the name for the file. The default is "biogeobears_data.txt".
#'
#' @return Writes a phylip format data file compatible with BioGeoBEARS.
#' @export
make_phylip_data <- function(data, states, file = "biogeobears_data.txt"){
    ## Function to produce phylip file to use in BioGeoBEARS.

    ## Check if file already exists and give a warning.
    if( file.exists(file.path(file)) ){
        stop( "File already exists. Will not overwrite it." )
    }

    outfile <- file(file.path(file), open="a")
    ## data is a named string with the states.

    ntips <- as.character( length( data ) )
    nstates <- as.character( length( states ) )
    states_chunk <- paste0( "(", paste(states, collapse = " "), ")")
    head_string <- paste(ntips, nstates, states_chunk, sep = "\t")
    cat(head_string, file = outfile, append = FALSE) ## First line.
    cat("\n", file = outfile, append = TRUE)

    data_mat <- makeDataTips(X = data, states = states)
    for( i in 1:nrow(data_mat) ){
        st_string <- paste(as.character(data_mat[i,]), collapse = "")
        spp_string <- paste(rownames(data_mat)[i], st_string, sep = "\t")
        cat(spp_string, file = outfile, append = TRUE)
        cat("\n", file = outfile, append = TRUE)
    }

    close(con = outfile)
}
