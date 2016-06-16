test_that("That calculating prior probs of interaction pattern params works", {
    skip_on_cran()

    # create an example distribution
    set.seed(12345)

    intercepts <- rnorm(4, mean = 0, sd = 2)
    latent_pos <- array(data = rnorm(n = 32, mean = 0, sd = 2), dim = c(4, 4, 2))
    coefficients <- matrix(rnorm(n = 16, mean = 0, sd = 2), nrow = 4, ncol = 4)

    # first lets try without covariates
    result <- test_internal_functions(
        Test_Prior_Pobability_Of_I_P_Params = TRUE,
        intercepts = intercepts,
        coefficients = coefficients,
        latent_positions = latent_pos,
        intercept_prior_mean = 0,
        intercept_prior_standard_deviation = 4,
        coefficient_prior_mean = 0,
        coefficient_prior_standard_deviation = 4,
        latent_position_prior_mean = 0,
        latent_position_prior_standard_deviation = 4,
        using_coefficients = FALSE)

    expect_that(result, equals(sum(dnorm(intercepts, 0, 4, TRUE)) + sum(dnorm(latent_pos, 0, 4, TRUE))))

    # make sure that using coefficients works as well
    result2 <- test_internal_functions(
        Test_Prior_Pobability_Of_I_P_Params = TRUE,
        intercepts = intercepts,
        coefficients = coefficients,
        latent_positions = latent_pos,
        intercept_prior_mean = 0,
        intercept_prior_standard_deviation = 4,
        coefficient_prior_mean = 0,
        coefficient_prior_standard_deviation = 4,
        latent_position_prior_mean = 0,
        latent_position_prior_standard_deviation = 4,
        using_coefficients = TRUE)

    expect_that(result2, equals(sum(dnorm(intercepts, 0, 4, TRUE)) +
                                  sum(dnorm(coefficients, 0, 4, TRUE)) +
                                  sum(dnorm(latent_pos, 0, 4, TRUE))))
})
