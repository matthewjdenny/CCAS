test_that("That LDA contribution works", {
    skip_on_cran()

    # create an example distribution
    set.seed(12345)

    num_topics <- 5
    num_terms <- 10
    # five topics
    alpha <- 1
    temp <- runif(num_topics)
    alpha_m <- alpha * temp / sum(temp)

    # ten unique words
    beta <- 2
    temp2 <- runif(num_terms)
    beta_n <- beta * temp2 / sum(temp2)

    word_type_topic_counts <- matrix(round(rpois(50, lambda = 10)),
                                           nrow = num_terms,
                                           ncol = num_topics)
    topic_token_counts <- colSums(word_type_topic_counts)

    topic <- 2
    current_word_type <- 2
    current_document_topic_counts <- c(5, 10, 2, 6, 1)
    tokens_in_document <- 24
    current_token_topic_assignment <- 2

    # first lets try without covariates
    result <- test_internal_functions(
        Test_LDA_Contribution = TRUE,
        tokens_in_document = tokens_in_document,
        current_token_topic_assignment = current_token_topic_assignment,
        current_document_topic_counts = current_document_topic_counts,
        word_type_topic_counts = word_type_topic_counts,
        topic_token_counts = topic_token_counts,
        topic = topic,
        current_word_type = current_word_type,
        alpha_m = alpha_m,
        beta_n = beta_n,
        beta = beta) # we did this so we do not have to keep re-calculating

    current_dtc <- current_document_topic_counts[topic + 1]
    wttc <- word_type_topic_counts[current_word_type + 1, topic + 1]
    tc <- topic_token_counts[topic + 1]

    if (topic == current_token_topic_assignment) {
      current_dtc <- current_dtc - 1
      wttc <- wttc - 1
      tc <- tc - 1
    }

    contribution <- log(current_dtc + alpha_m[topic + 1]) +
      log(wttc + beta_n[current_word_type + 1]) - log(tc + beta)

    expect_that(result, equals(contribution))
})
