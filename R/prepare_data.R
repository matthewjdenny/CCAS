#' @title Prepare communication network data for analysis.
#' @description Reads in and stores communication network data in a ComNet
#' object for use in main estimation procedure.
#'
#' @param document_authors A vector recording the index of the sender of each
#' document.
#' @param document_edge_matrix A documents x actors matrix recording whether
#' actor j was a recipient of message i.
#' @param document_term_matrix A documents x unique terms matrix recording the
#' count of unique term j in document i.
#' @param covariate_data Defaults to NULL.
#' @param vocabulary Defaults to NULL.
#' @return An object of class ComNet
prepare_data <- function(document_authors,
                         document_edge_matrix,
                         document_term_matrix,
                         covariate_data = NULL,
                         vocabulary = NULL){

    # initialize boolean indicating whether covariate data was provided
    using_covariates <- FALSE
    if (!is.null(covariate_data)) {
        using_covariates <- TRUE
    }


}
