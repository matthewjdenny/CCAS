initialize_LSM_priors <- function(CCAS_Object,
                                  LSM_intercept_prior_variance,
                                  LSM_intercept_prior_mean,
                                  LSM_position_prior_variance,
                                  LSM_position_prior_mean,
                                  LSM_coefficient_prior_variance,
                                  LSM_coefficient_prior_mean) {

    # currently just set them, could be different for different clusters in the
    # future
    CCAS_Object@LSM_intercept_prior_standard_deviation = sqrt(LSM_intercept_prior_variance)
    CCAS_Object@LSM_intercept_prior_mean = LSM_intercept_prior_mean
    CCAS_Object@LSM_position_prior_standard_deviation = sqrt(LSM_position_prior_variance)
    CCAS_Object@LSM_position_prior_mean = LSM_position_prior_mean
    CCAS_Object@LSM_coefficient_prior_standard_deviation = sqrt(LSM_coefficient_prior_variance)
    CCAS_Object@LSM_coefficient_prior_mean = LSM_coefficient_prior_mean

    return(CCAS_Object)
}

