initialize_latent_variables <- function(CCAS_Object){

    # first initialize instances of LSM variables centered at the prior
    # mean and variance
    using_coefficients <- CCAS_Object@ComNet_Object@using_covariates
    intercepts <- rep(CCAS_Object@LSM_intercept_prior_mean,
                      CCAS_Object@interaction_patterns)
    if (using_coefficients) {
        coefficients <- matrix(CCAS_Object@LSM_coefficient_prior_mean,
                               nrow = CCAS_Object@interaction_patterns,
                               ncol = CCAS_Object@number_of_covariates)
    } else {
        coefficients <- matrix(0,
                               nrow = CCAS_Object@interaction_patterns,
                               ncol = 2)
    }
    latent_positions <- array(CCAS_Object@LSM_position_prior_mean,
                              dim = c(CCAS_Object@interaction_patterns,
                                      CCAS_Object@ComNet_Object@num_actors,
                                      CCAS_Object@latent_space_dimensions))

    # need to make these into vectors since they are only a single value so far.
    # eventually could end up being different for each interaction pattern
    ipv <- rep(CCAS_Object@LSM_intercept_prior_standard_deviation,
               CCAS_Object@interaction_patterns)
    cpv <- rep(CCAS_Object@LSM_coefficient_prior_standard_deviation,
               CCAS_Object@interaction_patterns)
    ppv <- rep(CCAS_Object@LSM_position_prior_standard_deviation,
               CCAS_Object@interaction_patterns)

    # sample the new parameter values for these values centered at the prior
    # mean. We are using the mjd::sample_new_interaction_pattern_parameters
    # function here, which exposed to R through the snipp() wrapper.
    params1 <- snipp(
        intercepts,
        coefficients,
        latent_positions,
        intercept_proposal_standard_deviations = ipv,
        coefficient_proposal_standard_deviations = cpv,
        latent_position_proposal_standard_deviations = ppv,
        using_coefficients)

    # extract the initialized interaction pattern parameters from the returned
    # list and put them in a named list
    LSM_Params <- list(intercepts = params1[[1]],
                       coefficients = params1[[2]],
                       positions = params1[[3]])

    # now initialize topic interaction patterns (subtract one to be zero
    # indexed) as everything is zero indexed in c++
    topic_interaction_patterns <- sample(1:CCAS_Object@interaction_patterns,
                                         size = CCAS_Object@number_of_topics,
                                         replace = TRUE) - 1

    # now initialize token topic assignments using the
    # mjd::sample_token_topics_from_generative_process c++ function to
    # initialize then according to the generative process
    params <- sttgp(
        CCAS_Object@ComNet_Object@token_topic_assignment_list_zero_indexed,
        CCAS_Object@ComNet_Object@token_word_type_list_zero_indexed,
        CCAS_Object@alpha_m,
        CCAS_Object@beta_n,
        CCAS_Object@ComNet_Object@num_documents,
        FALSE,
        runif(n = 5*CCAS_Object@ComNet_Object@num_tokens))
    # for the last argument, just pass in lots of random numbers to be safe.

    # stick everything in a list object with the right names so it can be
    # returned
    LDA_Params <- list(token_topic_assignments = params[[1]],
                      token_word_types = params[[2]],
                      document_topic_counts = params[[3]],
                      topic_token_counts = params[[4]],
                      word_type_topic_counts = params[[5]],
                      document_topic_distributions = params[[6]],
                      topic_word_type_distributions = params[[7]],
                      topic_interaction_patterns_zero_indexed = topic_interaction_patterns)

    return(list(LSM_Params = LSM_Params,
                LDA_Params = LDA_Params))
}
