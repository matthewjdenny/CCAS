test_that("That we get it right", {
    skip_on_cran()
    skip() # until we get it working

    # we need a very simple model with a small number of actors, topics, odcuments etc.
    # create an example distribution
    set.seed(12345)

    GiR_samples = 10000000
    num_documents = 5
    words_per_doc = 10
    num_topics = 4
    num_terms = 10
    num_actors = 4
    num_ip = 2
    num_ld = 2
    num_covar = 1
    total_word_count = num_documents * words_per_doc
    # five topics
    alpha = 1
    temp <- runif(num_topics)
    alpha_m = alpha * temp/sum(temp)

    # ten unique words
    beta = 2
    temp2 <- runif(num_terms)
    beta_n = beta * temp2/sum(temp2)

    author_indexes <- floor(runif(n = num_documents, min = 0, max = num_actors - 0.00001))

    covars <- array(data = rnorm(n = num_actors*num_actors*num_covar,0,2),
                    dim = c(num_actors, num_actors, num_covar))



    # first lets try without covariates
    forward_samples <- test_internal_functions(
        Getting_It_Right = TRUE,
        author_indexes = author_indexes,
        covariates = covars,
        alpha_m = alpha_m,
        beta_n = beta_n,
        using_coefficients = TRUE,
        intercept_prior_mean = 0,
        intercept_prior_standard_deviation = 1,
        intercept_proposal_standard_deviations = c(0.25,0.25),
        coefficient_prior_mean = 0,
        coefficient_prior_standard_deviation = 1,
        coefficient_proposal_standard_deviations = c(0.25,0.25),
        latent_position_prior_mean = 0,
        latent_position_prior_standard_deviation =1,
        latent_position_proposal_standard_deviations = c(0.25,0.25),
        target_accept_rate = 0.25,
        tollerance = 0.05,
        update_size = 0.05,
        seed = 12345,
        iterations = 5,
        metropolis_iterations = 50,
        iterations_before_t_i_p_updates = 2,
        update_t_i_p_every_x_iterations = 2,
        perform_adaptive_metropolis = TRUE,
        slice_sample_alpha_m = -2,
        slice_sample_step_size = 1,
        parallel = FALSE,
        num_documents = num_documents,
        words_per_doc = words_per_doc,
        num_topics = num_topics,
        num_terms = num_terms,
        num_actors = num_actors,
        num_ip = num_ip,
        num_ld = num_ld,
        total_number_of_tokens = total_word_count,
        GiR_samples = GiR_samples,
        forward_sample = TRUE)

    backward_samples <- test_internal_functions(
        Getting_It_Right = TRUE,
        author_indexes = author_indexes,
        covariates = covars,
        alpha_m = alpha_m,
        beta_n = beta_n,
        using_coefficients = TRUE,
        intercept_prior_mean = 0,
        intercept_prior_standard_deviation = 1,
        intercept_proposal_standard_deviations = c(0.25,0.25),
        coefficient_prior_mean = 0,
        coefficient_prior_standard_deviation = 1,
        coefficient_proposal_standard_deviations = c(0.25,0.25),
        latent_position_prior_mean = 0,
        latent_position_prior_standard_deviation =1,
        latent_position_proposal_standard_deviations = c(0.25,0.25),
        target_accept_rate = 0.25,
        tollerance = 0.05,
        update_size = 0.05,
        seed = 12345,
        iterations = 5,
        metropolis_iterations = 50,
        iterations_before_t_i_p_updates = 2,
        update_t_i_p_every_x_iterations = 2,
        perform_adaptive_metropolis = TRUE,
        slice_sample_alpha_m = -2,
        slice_sample_step_size = 1,
        parallel = FALSE,
        num_documents = num_documents,
        words_per_doc = words_per_doc,
        num_topics = num_topics,
        num_terms = num_terms,
        num_actors = num_actors,
        num_ip = num_ip,
        num_ld = num_ld,
        total_number_of_tokens = total_word_count,
        GiR_samples = GiR_samples,
        forward_sample = FALSE)

    # now we need to compare

})
