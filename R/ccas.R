#' @title Infer Communication Content and Structure
#' @description Performs inference on the content conditional structure of a
#' text valued communication network.
#'
#' @param formula A formula object of the form 'ComNet ~
#' euclidean(d = 2)' where d is the number of dimensions in the latent space
#' that the user would like to include, and Comnet is a ComNet object generated
#' by the prepare_data() function. May also include optional terms
#' 'sender("covariate_name")',
#' 'receiver("covariate_name")' and 'nodemix("covariate_name", base = value)'
#' which are defined analogously to the arguments in the latentnet package.
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
ccas <- function(formula,
                 interaction_patterns = 4,
                 topics = 40,
                 alpha = 1,
                 beta = 0.01,
                 iterations = 1000,
                 metropolis_hastings_iterations = 500,
                 target_accept_rate = 0.25,
                 tollerance = 0.05) {

    # possible terms for inclusion in model specification.
    possible_structural_terms <- c("euclidean")
    possible_covariate_terms <- c("sender","receiver","nodemix")
    possible_network_terms <- c("netcov")

    # make sure that formula is a formula object
    formula <- as.formula(formula)

    # parse the forumla and return a list object containing the ComNet object
    # and the latent space model spcification
    parsed_specifcation <- parse_specification(formula,
                                              possible_structural_terms,
                                              possible_covariate_terms,
                                              possible_network_terms,
                                              terms_to_parse = "structural")

    # get the ComNet object from the specification
    ComNet_Object <- parsed_specifcation$ComNet

    # now extract indicator of whether we are using covariates
    using_covariates <- parsed_specifcation$using_covariates

    # if we are using covariates, then generate the covariate array:
    if (using_covariates) {
        temp <- generate_covariate_array(formula,
                                         possible_structural_terms,
                                         possible_covariate_terms,
                                         possible_network_terms,
                                         ComNet_Object)
        covariate_array <- temp$covariate_array
        number_of_covariates <- temp$number_of_covariates
    }


    # initialize an object of class CCAS to store everything. This will include
    # initializing all latent variables and organizing data in a format that is
    # appropriate for the main inference function



    # initialize all latent variable values
    # latent_variables <- initialize_latent_variables(CCAS_Object)

    # run inference

    # run MH to convergence

    # generate diagnostics

    # retrun the CCAS object
    return(NULL)
}
