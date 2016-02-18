ccas <- function(specification,
                 document_authors,
                 document_edge_matrix,
                 document_term_matrix,
                 covariate_data = NULL,
                 vocabulary = NULL,
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
    possible_terms <- c("euclidean","nodemix","sender","receiver","topic_model_parameters")

    # parse the specification and return a list object.
    parsed_specifcation <- parsed_specifcation(specification,
                                               using_covariates,
                                               covariate_data,
                                               possible_terms)

    # the following terms will be included in the specification string.
    # interaction_patterns = 4,
    # latent_space_dimensions = 2,
    # alpha = 1,
    # beta = 0.01

    # initialize an object of class CCAS to store everything. This will include
    # initializing all latent variables and organizing data in a format that is
    # appropriate for the main inference function

    # run inference

    # run MH to convergence

    # generate diagnostics

    # retrun the CCAS object
    return(NULL)
}
