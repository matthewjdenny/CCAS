test_that("That Adaptive_Metropolis works", {
    skip_on_cran()

    # create an example distribution
    set.seed(12345)

    accept_rates <- c(0.21,0.8,0.1,0.29)
    # first lets try without covariates
    result <- test_internal_functions(
        Test_Adaptive_Metropolis = TRUE,
        intercept_proposal_standard_deviations = c(0.5,0.5,0.5,0.5),
        coefficient_proposal_standard_deviations = c(0.5,0.5,0.5,0.5),
        latent_position_proposal_standard_deviations = c(0.5,0.5,0.5,0.5),
        accept_rates = accept_rates,
        target_accept_rate = 0.25,
        tollerance = 0.05,
        update_size = 0.05)

    res <- as.numeric(result[[1]])
    expect_equal(res, c(0.5,0.55,0.45,0.5))
    # no errors!

})
