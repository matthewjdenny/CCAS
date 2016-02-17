test_that("That calculating prior probs of interaction pattern params works", {
    skip_on_cran()

    # create an example distribution
    set.seed(12345)

    intercepts <- rnorm(4,mean = 0, sd = 2)
    latent_pos <- array(data = rnorm(n = 32, mean = 0, sd = 2), dim = c(4,4,2))
    coefficients <- matrix(rnorm(n = 16, mean = 0, sd = 2),nrow = 4, ncol = 4)

    # first lets try without covariates
    result <- test_internal_functions(
        Test_Prior_Pobability_Of_I_P_Params = TRUE,
        intercepts = intercepts,
        coefficients = coefficients,
        latent_pos = latent_pos,
        intercept_prior_mean = 0,
        intercept_prior_variance = 4,
        coefficient_prior_mean = 0,
        coefficient_prior_variance = 4,
        latent_position_prior_mean = 0,
        latent_position_prior_variance = 4,
        using_coefficients = FALSE)

    # no errors need to write an analytical test in R at some point

    # for now, I am just going to make sure it does not change from what I have
    # seen
    previous_result <- -87.7845433508862
    expect_equal(result,previous_result)


    # make sure that using coefficients works as well
    result2 <- test_internal_functions(
        Test_Prior_Pobability_Of_I_P_Params = TRUE,
        intercepts = intercepts,
        coefficients = coefficients,
        latent_pos = latent_pos,
        intercept_prior_mean = 0,
        intercept_prior_variance = 4,
        coefficient_prior_mean = 0,
        coefficient_prior_variance = 4,
        latent_position_prior_mean = 0,
        latent_position_prior_variance = 4,
        using_coefficients = TRUE)

})
