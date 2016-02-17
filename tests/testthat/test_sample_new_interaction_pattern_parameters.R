test_that("That sampling new interaction pattern params works", {
    skip_on_cran()

    # create an example distribution
    set.seed(12345)

    intercepts <- rnorm(4,mean = 0, sd = 2)
    latent_pos <- array(data = rnorm(n = 32, mean = 0, sd = 2), dim = c(4,4,2))
    coefficients <- matrix(rnorm(n = 16, mean = 0, sd = 2),nrow = 4, ncol = 4)

    # first lets try without covariates
    result <- test_internal_functions(
        Test_Sample_New_I_P_Parameters = TRUE,
        intercepts = intercepts,
        coefficients = coefficients,
        latent_pos = latent_pos,
        intercept_proposal_variances = c(0.5,0.5,0.5,0.5),
        coefficient_proposal_variances = c(0.5,0.5,0.5,0.5),
        latent_position_proposal_variances = c(0.5,0.5,0.5,0.5),
        using_coefficients = FALSE)

    # no errors, we should write a Monte Carlo test here against an R version

    # make sure that using coefficients works as well
    result2 <- test_internal_functions(
        Test_Sample_New_I_P_Parameters = TRUE,
        intercepts = intercepts,
        coefficients = coefficients,
        latent_pos = latent_pos,
        intercept_proposal_variances = c(0.5,0.5,0.5,0.5),
        coefficient_proposal_variances = c(0.5,0.5,0.5,0.5),
        latent_position_proposal_variances = c(0.5,0.5,0.5,0.5),
        using_coefficients = TRUE)

})
