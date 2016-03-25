test_that("main the parallel token topic distribution updates return the same thing", {
    skip_on_cran()

    skip("We are using threading here which will not work on travis")
    # create an example distribution
    set.seed(12345)
    data(ComNet_data)

    net_covariate <- matrix(runif(100),10,10)

    # specify a formula that we will use for testing.
    formula <- ComNet_data ~ euclidean(d = 2) +
        sender("age") +
        nodemix("gender", base = "male") +
        netcov("net_covariate")

    system.time({
        CCAS_Object <- ccas(formula,
                            interaction_patterns = 4,
                            topics = 100,
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
                            slice_sample_step_size = 1)
    })
    # user  system elapsed
    # 24.147   0.535  24.851

    system.time({
        CCAS_Object2 <- ccas(formula,
                             interaction_patterns = 4,
                             topics = 100,
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
                             parallel = TRUE)
    })
    # user  system elapsed
    # 25.094   0.548  19.262

    # extract the topic model results
    TM1 <- CCAS_Object@topic_model_results
    TM2 <- CCAS_Object2@topic_model_results

    # these should all be identical and all check out!
    expect_equal(TM1$topic_interaction_patterns,TM2$topic_interaction_patterns)
    expect_equal(TM1$topic_token_counts,TM2$topic_token_counts)
    expect_equal(TM1$word_type_topic_counts,TM2$word_type_topic_counts)
    expect_equal(TM1$token_topic_assignments,TM2$token_topic_assignments)
    expect_equal(TM1$document_topic_counts,TM2$document_topic_counts)


})
