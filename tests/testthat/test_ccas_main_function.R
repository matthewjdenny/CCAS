test_that("main function works", {
    skip_on_cran()

    # create an example distribution
    set.seed(12345)
    data(ComNet_data)

    net_covariate <- matrix(runif(100),10,10)

    # specify a formula that we will use for testing.
    formula <- ComNet_data ~ euclidean(d = 2) +
                             sender("age") +
                             nodemix("gender", base = "male") +
                             netcov("net_covariate")

    CCAS_Object <- ccas(formula,
                     interaction_patterns = 4,
                     topics = 40,
                     alpha = 1,
                     beta = 0.01,
                     iterations = 20,
                     metropolis_hastings_iterations = 500,
                     final_metropolis_hastings_burnin = 500,
                     target_accept_rate = 0.25,
                     tollerance = 0.05,
                     adaptive_metropolis_update_size = 0.05,
                     LSM_intercept_proposal_variance = .2,
                     LSM_position_proposal_variance = .2,
                     LSM_coefficient_proposal_variance = .2)

})
