test_that("That Update_Single_Token_Topic_Assignment works", {
    skip_on_cran()

    # create an example distribution
    set.seed(12345)

    edge_probs <- array(data = runif(n = 64), dim = c(4,4,4))

    num_topics = 5
    num_terms = 10
    # five topics
    alpha = 1
    temp <- runif(num_topics)
    alpha_m = alpha * temp/sum(temp)

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
        Test_Update_Single_Token_Topic_Assignment = TRUE,
        edge_probabilities = edge_probs,
        tokens_in_document = 24,
        current_token_topic_assignment = 2,
        current_document_topic_counts = c(5,10,2,6,1),
        word_type_topic_counts = word_type_topic_counts,
        topic_token_counts = topic_token_counts,
        topic = 2,
        current_word_type = 2,
        alpha_m = alpha_m,
        beta_n = beta_n,
        document_edge_values = c(0,0,1,1),
        topic_interaction_patterns = c(0,1,1,2,3),
        document_sender = 0 ,
        rand_num = runif(1))

    # no errors, will write an analytical test

})
