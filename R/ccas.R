#' @title Infer Communication Content and Structure
#' @description Performs inference on the content conditional structure of a
#' text valued communication network.
#'
#' @param specification A formulat object of the form 'ComNet ~
#' euclidean(d = 2)' where d is the number of dimensions in the latent space
#' that the user would like to include, and Comnet is a ComNet object generated
#' by the prepare_data() function. May also include optional terms
#' 'sender("covariate_name")',
#' 'receiver("covariate_name")' and 'nodemix("covariate_name", base = value)'
#' which are defined analogously to the arguments in the latentnet package.
#' @param document_authors A vector recording the index of the sender of each
#' document.
#' @param document_edge_matrix A documents x actors matrix recording whether
#' actor j was a recipient of message i.
#' @param document_term_matrix A documents x unique terms matrix recording the
#' count of unique term j in document i.
#' @param covariate_data Defaults to NULL.
#' @param vocabulary Defaults to NULL.
#' @param interaction_patterns Defaults to 4.
#' @param topics Defaults to 40.
#' @param alpha Defaults to 1.
#' @param beta Defaults to 0.01.
#' @param iterations Defaults to 1000.
#' @param metropolis_hastings_iterations Defaults to 500.
#' @param target_accept_rate Defaults to 0.25.
#' @param tollerance Defaults to 0.05.
#' @return An object of class CCAS containing estimation results.
#' @export
ccas <- function(specification,
                 document_authors,
                 document_edge_matrix,
                 document_term_matrix,
                 covariate_data = NULL,
                 vocabulary = NULL,
                 interaction_patterns = 4,
                 topics = 40,
                 alpha = 1,
                 beta = 0.01,
                 iterations = 1000,
                 metropolis_hastings_iterations = 500,
                 target_accept_rate = 0.25,
                 tollerance = 0.05) {

    # initialize boolean indicating whether covariate data was provided
    using_covariates <- FALSE
    if (!is.null(covariate_data)) {
        using_covariates <- TRUE
    }

    # possible terms for inclusion in model specification.
    possible_structural_terms <- c("euclidean")
    possible_covariate_terms <- c("sender","receiver","nodemix")
    possible_network_terms <- c("netcov")

    # make sure that specification is a formula object
    specification <- as.formula(specification)

    # parse the specification and return a list object.
    parsed_specifcation <- parsed_specifcation(specification,
                                               using_covariates,
                                               covariate_data,
                                               possible_structural_terms,
                                               possible_covariate_terms,
                                               possible_network_terms)

    # initialize an object of class CCAS to store everything. This will include
    # initializing all latent variables and organizing data in a format that is
    # appropriate for the main inference function

    # run inference

    # run MH to convergence

    # generate diagnostics

    # retrun the CCAS object
    return(NULL)
}
