test_that("main function works", {
    skip_on_cran()

    skip("Something weird is going on with Travis")
    # create an example distribution
    set.seed(12345)
    data(ComNet_data)

    # specify a formula that we will use for testing.
    formula <- ComNet_data ~ euclidean(d = 2) +
                             nodemix("Gender", base = "Male")

    CCAS_Object <- ccas(formula,
                     interaction_patterns = 4,
                     topics = 40,
                     alpha = 1,
                     beta = 0.01,
                     iterations = 20,
                     metropolis_hastings_iterations = 500,
                     final_metropolis_hastings_iterations = 10000,
                     final_metropolis_hastings_burnin = 5000,
                     thin = 1/10,
                     target_accept_rate = 0.25,
                     tolerance = 0.05,
                     adaptive_metropolis_update_size = 0.05,
                     LSM_proposal_variance = .5,
                     LSM_prior_variance = 1,
                     LSM_prior_mean = 0,
                     slice_sample_alpha_m = TRUE,
                     slice_sample_step_size = 1,
                     generate_plots = FALSE,
                     output_directory = NULL,
                     output_name_stem = NULL)

    expect_equal(nrow(CCAS_Object@MCMC_output$intercepts),1000)

})
