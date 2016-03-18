#' @title Infer Communication Content and Structure
#' @description Performs inference on the content conditional structure of a
#' text valued communication network.
#'
#' @param formula A formula object of the form 'ComNet ~
#' euclidean(d = 2)' where d is the number of dimensions in the latent space
#' that the user would like to include, and Comnet is a ComNet object generated
#' by the prepare_data() function. May also include optional terms
#' 'sender("covariate_name")', receiver("covariate_name")',
#' 'nodemix("covariate_name", base = value)' and netcov("network_covariate"),
#' which are defined analogously to the arguments in the latentnet package.
#' @param interaction_patterns Defaults to 4.
#' @param topics Defaults to 40.
#' @param alpha Defaults to 1.
#' @param beta Defaults to 0.01.
#' @param iterations Defaults to 1000.
#' @param metropolis_hastings_iterations Defaults to 500.
#' @param metropolis_hastings_burnin Defaults to 500.
#' @param target_accept_rate Defaults to 0.25.
#' @param tollerance Defaults to 0.05.
#' @param LSM_intercept_proposal_variance Defaults to 2.
#' @param LSM_intercept_prior_variance Defaults to 5.
#' @param LSM_intercept_prior_mean Defaults to 0.
#' @param LSM_position_proposal_variance Defaults to 2.
#' @param LSM_position_prior_variance Defaults to 5.
#' @param LSM_position_prior_mean Defaults to 0.
#' @param LSM_coefficient_proposal_variance Defaults to 2.
#' @param LSM_coefficient_prior_variance Defaults to 5.
#' @param LSM_coefficient_prior_mean Defaults to 0.
#' @param iterations_before_t_i_p_updates Defaults to 5.
#' @param update_t_i_p_every_x_iterations Defaults to 5.
#' @param adaptive_metropolis Defaults to TRUE.
#' @param adaptive_metropolis_update_size Defaults to 0.05.
#' @param seed Defaults to 12345.
#' @param final_metropolis_hastings_iterations Defaults to 100000.
#' @param adaptive_metropolis_every_x_iterations Defaults to 1000.
#' @param stop_adaptive_metropolis_after_x_updates Defualts to 50.
#' @return An object of class CCAS containing estimation results.
#' @export
ccas <- function(formula,
                 interaction_patterns = 4,
                 topics = 40,
                 alpha = 1,
                 beta = 0.01,
                 iterations = 1000,
                 metropolis_hastings_iterations = 500,
                 metropolis_hastings_burnin = 500,
                 target_accept_rate = 0.25,
                 tollerance = 0.05,
                 LSM_intercept_proposal_variance = 2,
                 LSM_intercept_prior_variance = 5,
                 LSM_intercept_prior_mean = 0,
                 LSM_position_proposal_variance = 2,
                 LSM_position_prior_variance = 5,
                 LSM_position_prior_mean = 0,
                 LSM_coefficient_proposal_variance = 2,
                 LSM_coefficient_prior_variance = 5,
                 LSM_coefficient_prior_mean = 0,
                 iterations_before_t_i_p_updates = 5,
                 update_t_i_p_every_x_iterations = 5,
                 adaptive_metropolis = TRUE,
                 adaptive_metropolis_update_size = 0.05,
                 seed = 12345,
                 final_metropolis_hastings_iterations = 100000,
                 adaptive_metropolis_every_x_iterations = 1000,
                 stop_adaptive_metropolis_after_x_updates = 50) {

    # set the seed
    set.seed(seed)

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
    number_of_covariates <- 0
    covariate_array <- array(0,dim = c(10,10,10))
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
    CCAS_Object <- new("CCAS",
       ComNet_Object = ComNet_Object,
       interaction_patterns = interaction_patterns,
       number_of_topics = topics,
       latent_space_dimensions = parsed_specifcation$d,
       target_accept_rate = target_accept_rate,
       tollerance = tollerance,
       formula = formula,
       alpha = alpha,
       beta = beta,
       iterations = iterations,
       metropolis_hastings_iterations = metropolis_hastings_iterations,
       metropolis_hastings_burnin = metropolis_hastings_burnin,
       number_of_covariates = number_of_covariates,
       iterations_before_t_i_p_updates = iterations_before_t_i_p_updates,
       update_t_i_p_every_x_iterations = update_t_i_p_every_x_iterations,
       perform_adaptive_metropolis = adaptive_metropolis,
       adaptive_metropolis_update_size = adaptive_metropolis_update_size,
       seed = seed,
       final_metropolis_hastings_iterations = final_metropolis_hastings_iterations)

    CCAS_Object@covariate_array <- covariate_array

    # initialie vectors of prior variances
    CCAS_Object <- initialize_LSM_priors(CCAS_Object,
                                      LSM_intercept_prior_variance,
                                      LSM_intercept_prior_mean,
                                      LSM_position_prior_variance,
                                      LSM_position_prior_mean,
                                      LSM_coefficient_prior_variance,
                                      LSM_coefficient_prior_mean)

    # initialize vectors of proposal variances
    CCAS_Object@LSM_intercept_proposal_variance <- rep(
        LSM_intercept_proposal_variance,
        interaction_patterns
    )
    CCAS_Object@LSM_position_proposal_variance <- rep(
        LSM_position_proposal_variance,
        interaction_patterns
    )
    CCAS_Object@LSM_coefficient_proposal_variance <- rep(
        LSM_coefficient_proposal_variance,
        interaction_patterns
    )

    CCAS_Object@alpha_m <- rep(CCAS_Object@alpha/CCAS_Object@number_of_topics,
                   CCAS_Object@number_of_topics)
    CCAS_Object@beta_n <- rep(CCAS_Object@beta/CCAS_Object@ComNet_Object@vocabulary_size,
                  CCAS_Object@ComNet_Object@vocabulary_size)


    # initialize all latent variable values
    CCAS_Object@latent_variables <- initialize_latent_variables(CCAS_Object)

    # run inference
    Inference_Results <- model_inference(
        CCAS_Object@ComNet_Object@document_authors_zero_indexed,
        CCAS_Object@ComNet_Object@document_edge_matrix,
        CCAS_Object@latent_variables$LDA_Params$document_topic_counts,
        CCAS_Object@latent_variables$LDA_Params$topic_interaction_patterns_zero_indexed,
        CCAS_Object@latent_variables$LDA_Params$word_type_topic_counts,
        CCAS_Object@latent_variables$LDA_Params$topic_token_counts,
        CCAS_Object@latent_variables$LDA_Params$token_topic_assignments,
        CCAS_Object@latent_variables$LDA_Params$token_word_types,
        CCAS_Object@latent_variables$LSM_Params$intercepts,
        CCAS_Object@latent_variables$LSM_Params$coefficients,
        CCAS_Object@latent_variables$LSM_Params$positions,
        CCAS_Object@covariate_array,
        CCAS_Object@alpha_m,
        CCAS_Object@beta_n,
        CCAS_Object@ComNet_Object@using_covariates,
        CCAS_Object@LSM_intercept_prior_mean,
        CCAS_Object@LSM_intercept_prior_variance,
        CCAS_Object@LSM_intercept_proposal_variance,
        CCAS_Object@LSM_coefficient_prior_mean,
        CCAS_Object@LSM_coefficient_prior_variance,
        CCAS_Object@LSM_coefficient_proposal_variance,
        CCAS_Object@LSM_position_prior_mean,
        CCAS_Object@LSM_position_prior_variance,
        CCAS_Object@LSM_position_proposal_variance,
        CCAS_Object@target_accept_rate,
        CCAS_Object@tollerance,
        CCAS_Object@adaptive_metropolis_update_size,
        CCAS_Object@seed,
        CCAS_Object@iterations,
        CCAS_Object@metropolis_hastings_iterations,
        CCAS_Object@ComNet_Object@num_tokens,
        CCAS_Object@iterations_before_t_i_p_updates,
        CCAS_Object@update_t_i_p_every_x_iterations,
        CCAS_Object@perform_adaptive_metropolis)

    # make sure all of the output
    MCMC_Results <- list(intercepts = Inference_Results[[2]],
                         coefficients = Inference_Results[[3]],
                         latent_positions = Inference_Results[[4]],
                         accept_rates = Inference_Results[[9]],
                         final_proposal_variances = Inference_Results[[10]])

    Topic_Model_Results <- list(topic_interaction_patterns = Inference_Results[[1]],
                         document_topic_counts = Inference_Results[[5]],
                         word_type_topic_counts = Inference_Results[[6]],
                         topic_token_counts = Inference_Results[[7]],
                         token_topic_assignments = Inference_Results[[8]])

    # assign them to the slot in our CCAS Object
    CCAS_Object@MCMC_output <- MCMC_Results
    CCAS_Object@topic_model_results <- Topic_Model_Results

    cat("\n\n###################################################\n\n")
    cat("Main Inference Complete: Running MH to Convergence...\n\n")
    cat("###################################################\n\n")

    # select the last iteraction of output from MH from the main inference
    # to seed MH to convergence
    ints <- CCAS_Object@MCMC_output$intercepts[
        nrow(CCAS_Object@MCMC_output$intercepts),]

    coefs <- CCAS_Object@MCMC_output$coefficients[,,
        dim(CCAS_Object@MCMC_output$coefficients)[3]]

    ld <- CCAS_Object@latent_space_dimensions
    ind <- dim(CCAS_Object@MCMC_output$latent_positions)[3] - ld

    lat_pos <- CCAS_Object@MCMC_output$latent_positions[,,ind:(ind + ld)]

    # run MH to convergence
    final_mh_results <- mh_to_convergence(
        CCAS_Object@ComNet_Object@document_authors_zero_indexed,
        CCAS_Object@ComNet_Object@document_edge_matrix,
        CCAS_Object@topic_model_results$document_topic_counts,
        CCAS_Object@topic_model_results$topic_interaction_patterns,
        ints,
        coefs,
        lat_pos,
        CCAS_Object@covariate_array,
        CCAS_Object@ComNet_Object@using_covariates,
        CCAS_Object@LSM_intercept_prior_mean,
        CCAS_Object@LSM_intercept_prior_variance,
        CCAS_Object@MCMC_output$final_proposal_variances,
        CCAS_Object@LSM_coefficient_prior_mean,
        CCAS_Object@LSM_coefficient_prior_variance,
        CCAS_Object@MCMC_output$final_proposal_variances,
        CCAS_Object@LSM_position_prior_mean,
        CCAS_Object@LSM_position_prior_variance,
        CCAS_Object@MCMC_output$final_proposal_variances,
        CCAS_Object@target_accept_rate,
        CCAS_Object@tollerance,
        CCAS_Object@adaptive_metropolis_update_size,
        CCAS_Object@seed,
        CCAS_Object@final_metropolis_hastings_iterations,
        adaptive_metropolis_every_x_iterations,
        stop_adaptive_metropolis_after_x_updates)

    # now slot in last updates
    CCAS_Object@MCMC_output$intercepts = final_mh_results[[1]]
    CCAS_Object@MCMC_output$coefficients = final_mh_results[[2]]
    CCAS_Object@MCMC_output$latent_positions = final_mh_results[[3]]
    CCAS_Object@MCMC_output$final_proposal_variances = final_mh_results[[4]]

    # generate diagnostics

    # retrun the CCAS object
    return(CCAS_Object)
}
