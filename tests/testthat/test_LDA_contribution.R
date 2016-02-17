test_that("That LDA contribution works", {
    skip_on_cran()

    # create an example distribution
    set.seed(12345)

    num_topics = 5
    num_terms = 10
    # five topics
    alpha = 1
    temp <- runif(num_topics)
    alpha_m = alpha * temp2/sum(temp2)

    # ten unique words
    beta = 2
    temp2 <- runif(num_terms)
    beta_n = beta * temp2/sum(temp2)

    word_type_topic_counts <- matrix(round(rpois(50,lambda = 10)),
                                           nrow = num_terms,
                                           ncol = num_topics)
    topic_token_counts <- colSums(word_type_topic_counts)

    # first lets try without covariates
    result <- test_internal_functions(
        Test_LDA_Contribution = TRUE,
        tokens_in_document = 24,
        current_token_topic_assignment = 2,
        current_document_topic_counts = c(5,10,2,6,1),
        word_type_topic_counts = word_type_topic_counts,
        topic_token_counts = topic_token_counts,
        topic = 2,
        current_word_type = 2,
        alpha_m = alpha_m,
        beta_n = beta_n,
        beta = beta) # we did this so we do not have to keep re-calculating

    # no errors, will write an analytical test

})
