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

    # need to make these into vectors since they are only a single value
    # when they are used elsewhere
    ipv <- rep(CCAS_Object@LSM_intercept_prior_variance,
               CCAS_Object@interaction_patterns)
    cpv <- rep(CCAS_Object@LSM_coefficient_prior_variance,
               CCAS_Object@interaction_patterns)
    ppv <- rep(CCAS_Object@LSM_position_prior_variance,
               CCAS_Object@interaction_patterns)

    # sample the new parameter values for these values
    LSM_Params <- snipp(
        intercepts,
        coefficients,
        latent_positions,
        intercept_proposal_variances = ipv,
        coefficient_proposal_variances = cpv,
        latent_position_proposal_variances = ppv,
        using_coefficients)

    # now initialize token topic assignments

    alpha_m <- rep(CCAS_Object@alpha/CCAS_Object@number_of_topics,
                   CCAS_Object@number_of_topics)
    beta_n <- rep(CCAS_Object@beta/CCAS_Object@ComNet_Object@vocabulary_size,
                  CCAS_Object@ComNet_Object@vocabulary_size)
    params <- sttgp(
        CCAS_Object@ComNet_Object@token_topic_assignment_list_zero_indexed,
        CCAS_Object@ComNet_Object@token_word_type_list_zero_indexed,
        alpha_m,
        beta_n,
        CCAS_Object@ComNet_Object@num_documents,
        FALSE,
        runif(n = 5*CCAS_Object@ComNet_Object@num_tokens))

    LDA_Params <- list(token_topic_assignments = params[[1]],
                      token_word_types = params[[2]],
                      document_topic_counts = params[[3]],
                      topic_token_counts = params[[4]],
                      word_type_topic_counts = params[[5]],
                      document_topic_distributions = params[[6]],
                      topic_word_type_distributions = params[[7]])

    return(list(LSM_Params = LSM_Params,
                LDA_Params = LDA_Params))
}
