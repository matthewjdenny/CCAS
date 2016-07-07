test_that("That Update_Single_Token_Topic_Assignment works", {
    skip_on_cran()

    # create an example distribution
    seed <- 12345
    set.seed(seed)

    edge_probs <- array(data = runif(n = 64), dim = c(4,4,4))

    num_topics <- 5
    num_terms <- 10
    # five topics
    alpha <- 1
    temp <- runif(num_topics)
    alpha_m <- alpha * temp/sum(temp)

    # ten unique words
    beta <- 2
    temp2 <- runif(num_terms)
    beta_n <- beta * temp2 / sum(temp2)

    word_type_topic_counts <- matrix(round(rpois(50, lambda = 10)),
                                     nrow = num_terms,
                                     ncol = num_topics)
    topic_token_counts <- colSums(word_type_topic_counts)
    tokens_in_document <- 24
    current_token_topic_assignment <- 2
    current_document_topic_counts <- c(5, 10, 2, 6, 1)
    topic <- 2
    current_word_type <- 2
    document_edge_values <- c(0, 0, 1, 1)
    topic_interaction_patterns <- c(0, 1, 1, 2, 3)
    document_sender <- 0
    rand_num <- runif(1)

    # first lets try without covariates
    result <- test_internal_functions(
        Test_Update_Single_Token_Topic_Assignment = TRUE,
        edge_probabilities = edge_probs,
        tokens_in_document = tokens_in_document,
        current_token_topic_assignment = current_token_topic_assignment,
        current_document_topic_counts = current_document_topic_counts,
        word_type_topic_counts = word_type_topic_counts,
        topic_token_counts = topic_token_counts,
        topic = topic,
        current_word_type = current_word_type,
        alpha_m = alpha_m,
        beta_n = beta_n,
        document_edge_values = document_edge_values,
        topic_interaction_patterns = topic_interaction_patterns,
        document_sender = document_sender,
        rand_num = rand_num,
        parallel = FALSE,
        use_cached_token_topic_distribution = FALSE,
        cached_token_topic_distribution = rep(-1, num_topics))

    unnormalized_log_distribution <- vector("numeric", length(alpha_m))
    for (i in 1:length(alpha_m)) {
      lsm_contribution <- test_internal_functions(Test_LSM_Contribution = TRUE,
                                                  edge_probabilities = edge_probs,
                                                  tokens_in_document = tokens_in_document,
                                                  topic = i - 1,
                                                  current_token_topic_assignment = current_token_topic_assignment,
                                                  current_document_topic_counts = current_document_topic_counts,
                                                  document_edge_values = document_edge_values,
                                                  topic_interaction_patterns = topic_interaction_patterns,
                                                  document_sender = document_sender)
      lda_contribution <- test_internal_functions(Test_LDA_Contribution = TRUE,
                                                  tokens_in_document = tokens_in_document,
                                                  current_token_topic_assignment = current_token_topic_assignment,
                                                  current_document_topic_counts = current_document_topic_counts,
                                                  word_type_topic_counts = word_type_topic_counts,
                                                  topic_token_counts = topic_token_counts,
                                                  topic = i - 1,
                                                  current_word_type = current_word_type,
                                                  alpha_m = alpha_m,
                                                  beta_n = beta_n,
                                                  beta = sum(beta_n))
      unnormalized_log_distribution[i] <- lsm_contribution + lda_contribution
    }
    
    assignment <- test_internal_functions(
      Test_Log_Space_Multinomial_Sampler = TRUE,
      distribution = unnormalized_log_distribution, seed = seed, u = rand_num)

    result_parallel <- test_internal_functions(
        Test_Update_Single_Token_Topic_Assignment = TRUE,
        edge_probabilities = edge_probs,
        tokens_in_document = tokens_in_document,
        current_token_topic_assignment = current_token_topic_assignment,
        current_document_topic_counts = current_document_topic_counts,
        word_type_topic_counts = word_type_topic_counts,
        topic_token_counts = topic_token_counts,
        topic = topic,
        current_word_type = current_word_type,
        alpha_m = alpha_m,
        beta_n = beta_n,
        document_edge_values = document_edge_values,
        topic_interaction_patterns = topic_interaction_patterns,
        document_sender = document_sender,
        rand_num = rand_num,
        parallel = TRUE,
        use_cached_token_topic_distribution = FALSE,
        cached_token_topic_distribution = rep(-1, num_topics))


    expect_that(assignment - 1, equals(result[1]))
    expect_that(unnormalized_log_distribution, equals(result[-1]))

    expect_that(assignment - 1, equals(result_parallel[1]))
    expect_that(unnormalized_log_distribution, equals(result_parallel[-1]))
})
