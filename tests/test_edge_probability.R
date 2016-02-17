test_that("That calculating edge probability works", {
    skip_on_cran()

    # create an example distribution
    set.seed(12345)

    latent_pos <- array(data = rnorm(n = 18, mean = 0, sd = 10), dim = c(3,3,2))
    coefficients <- matrix(rnorm(n = 16, mean = 0, sd = 2),nrow = 4, ncol = 4)

    # first lets try without covariates
    result <- test_internal_functions(
        Test_Edge_Probability = TRUE,
        intercepts =  c(1,2,3,4),
        coefficients = coefficients,
        latent_pos = latent_pos,
        document_sender = 0, # note that these have to be zero indexed
        document_recipient = 2, # note that these have to be zero indexed
        current_covariates = c(100,100,100,100),
        interaction_pattern = 2, # note that these have to be zero indexed
        using_coefficients = FALSE)

})
