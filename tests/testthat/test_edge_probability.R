test_that("That calculating edge probability works", {
    skip_on_cran()

    # create an example distribution
    set.seed(12345)

    latent_pos <- array(data = rnorm(n = 32, mean = 0, sd = 2), dim = c(4,4,2))
    coefficients <- matrix(rnorm(n = 16, mean = 0, sd = 2),nrow = 4, ncol = 4)

    # first lets try without covariates
    result <- test_internal_functions(
        Test_Edge_Probability = TRUE,
        intercepts =  c(1.0,2.0,3.0,4.0),
        coefficients = coefficients,
        latent_pos = latent_pos,
        document_sender = 0, # note that these have to be zero indexed
        document_recipient = 2, # note that these have to be zero indexed
        current_covariates = c(0,2.0,0,-4.0),
        interaction_pattern = 2, # note that these have to be zero indexed
        using_coefficients = FALSE)

    # calculate the linear predictor by hand,
    # interaction pattern 3, actors 1 and 3
    lin_pred <- 3 - sqrt((latent_pos[3,1,1] - latent_pos[3,3,1])^2 +
                             (latent_pos[3,1,2] - latent_pos[3,3,2])^2)
    hand_calculation <- plogis(lin_pred)

    # I have noticed some small differences in the values of the hand and
    # function calculated values after the 4th decimal place. I have confirmed
    # that these result from rounding that gets done when passing numbers into
    # C++ using the NumericVector class in Rcpp. Will look to see if there is
    # another way we can pass stuff in to get exact results. Does not affect
    # functions that do not need to pass in arrays.

    expect_equal(round(result,4),round(hand_calculation,4))

    #now try with coefficients
    result2 <- test_internal_functions(
        Test_Edge_Probability = TRUE,
        intercepts =  c(1,2,3,4),
        coefficients = coefficients,
        latent_pos = latent_pos,
        document_sender = 0, # note that these have to be zero indexed
        document_recipient = 2, # note that these have to be zero indexed
        current_covariates = c(0,3,0,4),
        interaction_pattern = 2, # note that these have to be zero indexed
        using_coefficients = TRUE)

    # calculate the linear predictor by hand,
    # interaction pattern 3, actors 1 and 3
    lin_pred2 <- 3  + coefficients[3,2]*3 + coefficients[3,4]*4 -
        sqrt((latent_pos[3,1,1] - latent_pos[3,3,1])^2 +
        (latent_pos[3,1,2] - latent_pos[3,3,2])^2)
    hand_calculation2 <- plogis(lin_pred2)

    expect_equal(round(result2,4),round(hand_calculation2,4))

})
