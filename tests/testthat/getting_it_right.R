test_that("That we get it right", {
    skip_on_cran()

    # we need a very simple model with a small number of actors, topics, odcuments etc.
    # create an example distribution
    set.seed(12345)

    resample_token_word_types = FALSE
    GiR_samples = 5000000
    num_documents = 5
    words_per_doc = 4
    num_topics = 3
    num_terms = 5
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

    token_topic_assignments <- vector(mode = "list",
      length = num_documents)
    token_word_types <- vector(mode = "list",
                               length = num_documents)

    for (i in 1:num_documents) {
        # initialize to -1 to make sure we throw an error if the C++ generative
        # process does not take care of it.
        token_topic_assignments[[i]] <- rep(-1, words_per_doc)
        if (resample_token_word_types) {
            token_word_types[[i]] <- rep(-1, words_per_doc)
        } else {
            # if we are not resampling then they need to be fixed from the start.
            token_word_types[[i]] <- floor(runif(n = words_per_doc,
                                                 min = 0,
                                                 max = num_terms - .000001))
        }
    }


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
        forward_sample = TRUE,
        token_topic_assignments = token_topic_assignments,
        token_word_types = token_word_types,
        resample_word_types = resample_token_word_types,
        verbose = FALSE)

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
        forward_sample = FALSE,
        token_topic_assignments = token_topic_assignments,
        token_word_types = token_word_types,
        resample_word_types = resample_token_word_types,
        verbose = FALSE)

    # now we need to compare the two output streams
    # subsample points in the middle to make this possible to
    # plot for a really large number of samples
    k <- floor(0.0001 * GiR_samples)
    m <- .001 * GiR_samples
    e <- floor(0.0005 * GiR_samples)

    quant.subsample = function(y, m = 100, e = 1) {
      # m: size of a systematic sample
      # e: number of extreme values at either end to use
      x <- sort(y)
      n <- length(x)
      quants <- (1 + sin(1:m / (m + 1) * pi - pi / 2)) / 2
      sort(c(x[1:e], quantile(x, probs = quants), x[(n + 1 - e):n]))
      # Returns m + 2*e sorted values from the EDF of y
    }

    plt = do.call("rbind", lapply(seq_len(ncol(forward_samples)), function(x) {
      qq = data.frame("forward" = quant.subsample(forward_samples[, x], m, e),
        "backward" = quant.subsample(backward_samples[, x], m, e))
      ## qq = as.data.frame(qqplot(forward_samples[, x], backward_samples[, x],
      ##   plot.it = FALSE))
      qq$variable = colnames(forward_samples)[x]
      qq
    }))

    # we do not want the plot to auto-generate under travis and r cmd check
    make_plot <- FALSE
    if (make_plot) {
        pdf(file = "~/Desktop/QQ_Plots.pdf", height = 16, width = 20)
        ggplot2::ggplot(plt, ggplot2::aes(forward, backward)) + ggplot2::geom_point() +
            ggplot2::facet_wrap(~ variable, scales = "free")
        dev.off()
    }

})
