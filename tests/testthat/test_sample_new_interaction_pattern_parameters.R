test_that("That sampling new interaction pattern params works", {
    skip_on_cran()

    # create an example distribution
    seed <- 12345
    set.seed(seed)

    intercepts <- rnorm(1, mean = 0, sd = 2)
    latent_pos <- array(data = rnorm(n = 1, mean = 0, sd = 2), dim = c(1, 1, 1))
    coefficients <- matrix(rnorm(n = 1, mean = 0, sd = 2), nrow = 1, ncol = 1)

    set.seed(seed)
    # first lets try without covariates
    result <- test_internal_functions(
        Test_Sample_New_I_P_Parameters = TRUE,
        intercepts = intercepts,
        coefficients = coefficients,
        latent_positions = latent_pos,
        intercept_proposal_standard_deviations = .5,
        coefficient_proposal_standard_deviations = .5,
        latent_position_proposal_standard_deviations = .5,
        using_coefficients = FALSE)

    set.seed(seed)
    comp <- list(
      rnorm(1, intercepts, .5),
      rnorm(1, latent_pos, .5)
    )

    expect_that(as.numeric(result[[1]]), equals(comp[[1]]))
    expect_that(as.numeric(result[[3]]), equals(comp[[2]]))

    set.seed(seed)
    # make sure that using coefficients works as well
    result2 <- test_internal_functions(
        Test_Sample_New_I_P_Parameters = TRUE,
        intercepts = intercepts,
        coefficients = coefficients,
        latent_positions = latent_pos,
        intercept_proposal_standard_deviations = .5,
        coefficient_proposal_standard_deviations = .5,
        latent_position_proposal_standard_deviations = .5,
        using_coefficients = TRUE)

    set.seed(seed)
    comp <- list(
      rnorm(1, intercepts, .5),
      rnorm(1, coefficients, .5),
      rnorm(1, latent_pos, .5)
    )

    expect_that(as.numeric(result2[[1]]), equals(comp[[1]]))
    expect_that(as.numeric(result2[[2]]), equals(comp[[2]]))
    expect_that(as.numeric(result2[[3]]), equals(comp[[3]]))
})
