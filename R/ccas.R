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
#' @param interaction_patterns The number of different interaction patterns
#' governing message sending and recieving under the model. Defaults to 4.
#' @param topics The number of topics to be used in the model. Defaults to 40.
#' @param alpha The hyperparameter govering document-topic distributions. Lower
#' values encourage more peaked distributions. Defaults to 1.
#' @param beta The hyperparameter governing the Dirichlet prior on the
#' topic-word distributions. Lower values encourage more peaked distributions.
#' Defaults to 0.01.
#' @param iterations The number of iterations of Metropolis-within-Gibbs
#' sampling to be used in model estimation. Defaults to 1,000.
#' @param metropolis_hastings_iterations The number of Metropolis Hastings
#' iterations to be run durring each iteration of Metropolis-within-Gibbs
#' sampling to update interaction pattern parameters. Defaults to 500.
#' @param final_metropolis_hastings_burnin The number of iterations to of
#' Metropolis Hastings run after completing all main iterations of Gibbs
#' sampling to discard before keeping samples. Defaults to 50,000.
#' @param final_metropolis_hastings_iterations The number of iterations to run
#' Metropolis Hastings after completing all main iterations of Gibbs sampling.
#' Defaults to 100,000.
#' @param thin The proportion of network samples to keep from the final run of
#' Metropolis Hastings to convergence. Defaults to 1/100, meaning that every
#' 100'th network sample will be returned.
#' @param target_accept_rate The target acceptance rate for the Metropolis
#' Hastings algorithm. Defaults to 0.25 which is standard in the literature.
#' @param tolerance The tolerance for differences between the observed and
#' target Metropolis Hastings accept rates (+-). Defaults to 0.05.
#' @param LSM_proposal_variance The Metropolis Hastings proposal variance for
#' all interaction pattern parameters. Defaults to .5.
#' @param LSM_prior_variance The variance of the multivariate normal prior on
#' all interaction pattern parameters. Defaults to 1.
#' @param LSM_prior_mean The mean of the multivariate normal prior on all
#' interaction pattern parameters. Defaults to 0.
#' @param iterations_before_t_i_p_updates The number of iterations to wait
#' before beginning updates to topic interaction pattern assignments. Defaults
#' to 5.
#' @param update_t_i_p_every_x_iterations The number of iterations between
#' updates to topic interaction pattern assignments. Defaults to 5.
#' @param adaptive_metropolis Logical indicating whether adaptive Metropolis
#' should be used (whether the proposal variance should be optimized). Defaults
#' to TRUE.
#' @param adaptive_metropolis_update_size The amount by which the MH proposal
#' variance is changed (up or down) durring adaptive Metropolis. Defaults to
#' 0.05.
#' @param seed The seed to be used (for replicability across runs). Defaults to
#' 12345.
#' @param adaptive_metropolis_every_x_iterations The nubmer of iterations
#' between proposal variance updates durring the final run of MH to convergence.
#' Defaults to 1000.
#' @param stop_adaptive_metropolis_after_x_updates The number of Metropolis
#' Hastings proposal variance updates to complete durring the final run of
#' MH to convergence before fixing its value. Defualts to 50.
#' @param slice_sample_alpha_m Logical indicating whether hyperparameter
#' optimization should be used to determine the optimal value of alpha. Defaults
#' to FALSE.
#' @param slice_sample_step_size The initial size of the slice to use when slice
#' sampling alpha (hyperparameter optimization). Defaults to 1.
#' @param parallel Argument indicating whether the token topic distributions
#' should be generated in parallel. Defaults to FALSE. Can significantly reduce
#' runtime when training a model with a large number of topics.
#' @param cores The number of cores to be used if the parallel option is set to
#' TRUE. Should not exceed the nubmer of cores availabel on the machine and will
#' not show performance gains if cores > topics.
#' @return An object of class CCAS containing estimation results.
#' @export
ccas <- function(formula,
                 interaction_patterns = 4,
                 topics = 40,
                 alpha = 1,
                 beta = 0.01,
                 iterations = 1000,
                 metropolis_hastings_iterations = 500,
                 final_metropolis_hastings_burnin = 50000,
                 final_metropolis_hastings_iterations = 100000,
                 thin = 1/100,
                 target_accept_rate = 0.25,
                 tolerance = 0.05,
                 LSM_proposal_variance = .5,
                 LSM_prior_variance = 1,
                 LSM_prior_mean = 0,
                 iterations_before_t_i_p_updates = 5,
                 update_t_i_p_every_x_iterations = 5,
                 adaptive_metropolis = TRUE,
                 adaptive_metropolis_update_size = 0.05,
                 seed = 12345,
                 adaptive_metropolis_every_x_iterations = 1000,
                 stop_adaptive_metropolis_after_x_updates = 50,
                 slice_sample_alpha_m = FALSE,
                 slice_sample_step_size = 1,
                 parallel = FALSE,
                 cores = 2) {

    # for now, we only allow a common proposal variance, prior mean, and prior
    # variance. In the future, these could be different for each cluster and
    # parameter type.
    LSM_intercept_proposal_variance = LSM_proposal_variance
    LSM_intercept_prior_variance = LSM_prior_variance
    LSM_intercept_prior_mean = LSM_prior_mean
    LSM_position_proposal_variance = LSM_proposal_variance
    LSM_position_prior_variance = LSM_prior_variance
    LSM_position_prior_mean = LSM_prior_mean
    LSM_coefficient_proposal_variance = LSM_proposal_variance
    LSM_coefficient_prior_variance = LSM_prior_variance
    LSM_coefficient_prior_mean = LSM_prior_mean

    # set the seed
    set.seed(seed)

    # set the number of threads to use with parallel
    if (parallel) {
        RcppParallel::setThreadOptions(numThreads = cores)
    }

    # possible terms for inclusion in model specification.
    possible_structural_terms <- c("euclidean")
    possible_covariate_terms <- c("sender","receiver","nodemix")
    possible_network_terms <- c("netcov")

    # make sure that formula is a formula object
    formula <- as.formula(formula)

    # if we are slice sampling, then set the inteval to 5, otherwise, it will
    # be negative 2
    slice_sample_alpha_m_every <- -2
    if (slice_sample_alpha_m) {
        slice_sample_alpha_m_every <- 5
    }

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
       tolerance = tolerance,
       formula = formula,
       alpha = alpha,
       beta = beta,
       iterations = iterations,
       metropolis_hastings_iterations = metropolis_hastings_iterations,
       final_metropolis_hastings_burnin = final_metropolis_hastings_burnin,
       final_metropolis_hastings_iterations = final_metropolis_hastings_iterations,
       thin = thin,
       number_of_covariates = number_of_covariates,
       iterations_before_t_i_p_updates = iterations_before_t_i_p_updates,
       update_t_i_p_every_x_iterations = update_t_i_p_every_x_iterations,
       perform_adaptive_metropolis = adaptive_metropolis,
       adaptive_metropolis_update_size = adaptive_metropolis_update_size,
       seed = seed)

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
    CCAS_Object@LSM_intercept_proposal_standard_deviation <- rep(
        sqrt(LSM_intercept_proposal_variance),
        interaction_patterns
    )
    CCAS_Object@LSM_position_proposal_standard_deviation <- rep(
        sqrt(LSM_position_proposal_variance),
        interaction_patterns
    )
    CCAS_Object@LSM_coefficient_proposal_standard_deviation <- rep(
        sqrt(LSM_coefficient_proposal_variance),
        interaction_patterns
    )

    CCAS_Object@alpha_m <- rep(CCAS_Object@alpha/CCAS_Object@number_of_topics,
                   CCAS_Object@number_of_topics)
    CCAS_Object@beta_n <- rep(CCAS_Object@beta/CCAS_Object@ComNet_Object@vocabulary_size,
                  CCAS_Object@ComNet_Object@vocabulary_size)


    # initialize all latent variable values
    CCAS_Object@latent_variables <- initialize_latent_variables(CCAS_Object)

    cat("\n#############################################\n\n")
    cat("Initialization Complete: Running Inference...\n\n")
    cat("#############################################\n\n")

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
        CCAS_Object@LSM_intercept_prior_standard_deviation,
        CCAS_Object@LSM_intercept_proposal_standard_deviation,
        CCAS_Object@LSM_coefficient_prior_mean,
        CCAS_Object@LSM_coefficient_prior_standard_deviation,
        CCAS_Object@LSM_coefficient_proposal_standard_deviation,
        CCAS_Object@LSM_position_prior_mean,
        CCAS_Object@LSM_position_prior_standard_deviation,
        CCAS_Object@LSM_position_proposal_standard_deviation,
        CCAS_Object@target_accept_rate,
        CCAS_Object@tolerance,
        CCAS_Object@adaptive_metropolis_update_size,
        CCAS_Object@seed,
        CCAS_Object@iterations,
        CCAS_Object@metropolis_hastings_iterations,
        CCAS_Object@ComNet_Object@num_tokens,
        CCAS_Object@iterations_before_t_i_p_updates,
        CCAS_Object@update_t_i_p_every_x_iterations,
        CCAS_Object@perform_adaptive_metropolis,
        slice_sample_alpha_m_every,
        slice_sample_step_size,
        parallel)

    # make sure all of the output
    MCMC_Results <- list(intercepts = Inference_Results[[2]],
                         coefficients = Inference_Results[[3]],
                         latent_positions = Inference_Results[[4]],
                         accept_rates = Inference_Results[[9]],
                         final_proposal_standard_deviations = Inference_Results[[10]])

    Topic_Model_Results <- list(topic_interaction_patterns = Inference_Results[[1]],
                         document_topic_counts = Inference_Results[[5]],
                         word_type_topic_counts = Inference_Results[[6]],
                         topic_token_counts = Inference_Results[[7]],
                         token_topic_assignments = Inference_Results[[8]],
                         unnoramlized_LDA_log_likelihood = Inference_Results[[11]])

    # assign them to the slot in our CCAS Object
    CCAS_Object@MCMC_output <- MCMC_Results
    CCAS_Object@topic_model_results <- Topic_Model_Results

    cat("\n####################################################\n\n")
    cat("Main Inference Complete: Running MH to Convergence...\n\n")
    cat("####################################################\n\n")

    # select the last iteraction of output from MH from the main inference
    # to seed MH to convergence
    ints <- CCAS_Object@MCMC_output$intercepts[
        nrow(CCAS_Object@MCMC_output$intercepts),]
    coefs <- CCAS_Object@MCMC_output$coefficients[,,
        dim(CCAS_Object@MCMC_output$coefficients)[3]]
    ld <- CCAS_Object@latent_space_dimensions
    ind <- dim(CCAS_Object@MCMC_output$latent_positions)[3] - ld
    lat_pos <- CCAS_Object@MCMC_output$latent_positions[,,ind:(ind + ld)]

    # we also have to determine how many observations to actually store after
    # the burnin is complete and the chain has been appropriately thinned. We
    # do this outside of C++ since the C++ ceiling function can produce warnings
    # that will not pass R CMD Check.
    samples_to_store <- ceiling(CCAS_Object@thin *
        CCAS_Object@final_metropolis_hastings_iterations)

    # now calculate the inverse of the thin parameter to get our sample every
    # parameter
    sample_every <- (1/CCAS_Object@thin)

    if (sample_every < 1) {
        sample_every <- 1
        cat("You specified a value for thin > 1. You must select a value for thin < 1. All iterations will be kept...\n")
    }

    # the total iterations is the burning plus the iterations
    iters <- CCAS_Object@final_metropolis_hastings_iterations +
        CCAS_Object@final_metropolis_hastings_burnin

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
        CCAS_Object@LSM_intercept_prior_standard_deviation,
        CCAS_Object@MCMC_output$final_proposal_standard_deviations,
        CCAS_Object@LSM_coefficient_prior_mean,
        CCAS_Object@LSM_coefficient_prior_standard_deviation,
        CCAS_Object@MCMC_output$final_proposal_standard_deviations,
        CCAS_Object@LSM_position_prior_mean,
        CCAS_Object@LSM_position_prior_standard_deviation,
        CCAS_Object@MCMC_output$final_proposal_standard_deviations,
        CCAS_Object@target_accept_rate,
        CCAS_Object@tolerance,
        CCAS_Object@adaptive_metropolis_update_size,
        CCAS_Object@seed,
        iters,
        adaptive_metropolis_every_x_iterations,
        stop_adaptive_metropolis_after_x_updates,
        samples_to_store,
        sample_every,
        CCAS_Object@final_metropolis_hastings_burnin - 1
        )

    # now slot in last updates
    CCAS_Object@MCMC_output$intercepts = final_mh_results[[1]]
    CCAS_Object@MCMC_output$coefficients = final_mh_results[[2]]
    CCAS_Object@MCMC_output$latent_positions = final_mh_results[[3]]
    CCAS_Object@MCMC_output$final_proposal_standard_deviations = final_mh_results[[4]]

    cat("\n##############################\n\n")
    cat("Generating Diagnostic Plots...\n\n")
    cat("##############################\n\n")

    # generate diagnostics


    # retrun the CCAS object
    return(CCAS_Object)
}
