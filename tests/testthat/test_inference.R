test_that("That Inference works", {
    skip_on_cran()

    # create an example distribution
    set.seed(12345)

    num_documents = 20
    words_per_doc = 10
    num_topics = 5
    num_terms = 10
    num_actors = 4
    num_ip = 4
    numld = 2
    num_covar = 4
    total_word_count = num_documents * words_per_doc
    # five topics
    alpha = 1
    temp <- runif(num_topics)
    alpha_m = alpha * temp/sum(temp)

    # ten unique words
    beta = 2
    temp2 <- runif(num_terms)
    beta_n = beta * temp2/sum(temp2)

    intercepts <- rnorm(num_ip,mean = 0, sd = 2)
    latent_pos <- array(data = rnorm(n = num_ip*num_actors*numld, mean = 0, sd = 2),
                        dim = c(num_ip,num_actors,numld))
    coefficients <- matrix(rnorm(n = num_covar*num_ip, mean = 0, sd = 2),
                           nrow = num_covar, ncol = num_ip)

    document_edge_matrix <- matrix(round(runif(num_actors * num_documents)),
                                   nrow = num_documents,
                                   ncol = num_actors)

    author_indexes <- floor(runif(n = num_documents, min = 0, max = num_actors - 0.00001))

    covars <- array(data = rnorm(n = num_actors*num_actors*num_covar,0,2),
                    dim = c(num_actors, num_actors, num_covar))

    # create token topic assignemnts and token word types list
    init_token_topic_assignments <- vector(mode = "list",
                                           length = num_documents)
    token_word_types <- vector(mode = "list",
                               length = num_documents)
    document_topic_counts <- matrix(0,nrow = num_documents,
                                    ncol = num_topics)
    word_type_topic_counts <- matrix(0,
                                     nrow = num_terms,
                                     ncol = num_topics)
    for (i in 1:num_documents) {
        init_token_topic_assignments[[i]] <- floor(runif(n = words_per_doc,
                                                         min = 0,
                                                         max = num_topics - .000001))
        token_word_types[[i]] <- floor(runif(n = words_per_doc,
                                             min = 0,
                                             max = num_terms - .000001))
        for(j in 0:(num_topics-1)){
            document_topic_counts[i,j] <- length(which(init_token_topic_assignments[[i]]== j))
        }
        for(l in 1:words_per_doc){
            for(j in 1:num_terms){
                for(k in 1:num_topics){
                    if(token_word_types[[i]][l] == (j - 1) & init_token_topic_assignments[[i]][l] == (k -1)){
                        word_type_topic_counts[j,k] <- word_type_topic_counts[j,k] + 1
                    }
                }
            }
        }
    }

    topic_token_counts <- colSums(word_type_topic_counts)

    edge_probs <- array(data = runif(n = 64), dim = c(4,4,4))

    topic_interaction_patterns = c(0,1,1,2,3)

    # first lets try without covariates
    result <- test_internal_functions(
        Test_Inference = TRUE,
        author_indexes = author_indexes,
        document_edge_matrix = document_edge_matrix,
        document_topic_counts = document_topic_counts,
        topic_interaction_patterns = topic_interaction_patterns,
        word_type_topic_counts = word_type_topic_counts,
        topic_token_counts = topic_token_counts,
        token_topic_assignments = init_token_topic_assignments,
        token_word_types = token_word_types,
        intercepts = intercepts,
        coefficients = coefficients,
        latent_positions = latent_pos,
        covariates = covars,
        alpha_m = alpha_m,
        beta_n = beta_n,
        using_coefficients = TRUE,
        intercept_prior_mean = 0,
        intercept_prior_variance = 4,
        intercept_proposal_variances = c(0.5,0.5,0.5,0.5),
        coefficient_prior_mean = 0,
        coefficient_prior_variance = 4,
        coefficient_proposal_variances = c(0.5,0.5,0.5,0.5),
        latent_position_prior_mean = 0,
        latent_position_prior_variance= 4,
        latent_position_proposal_variances = c(0.5,0.5,0.5,0.5),
        target_accept_rate = 0.25,
        tollerance = 0.05,
        update_size = 0.05,
        seed = 12345,
        iterations = 20,
        metropolis_iterations = 10,
        total_number_of_tokens = total_word_count ,
        iterations_before_t_i_p_updates = 2,
        update_t_i_p_every_x_iterations = 2,
        perform_adaptive_metropolis = TRUE)

    # no errors!

})
