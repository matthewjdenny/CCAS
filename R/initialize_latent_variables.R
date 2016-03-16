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


    return(list(LSM_Params = LSM_Params,
                LDA_Params = NULL))
}
