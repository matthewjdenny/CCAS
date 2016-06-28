test_that("That LSM contribution works", {
    skip_on_cran()

    # create an example distribution
    set.seed(12345)

    edge_probs <- array(data = runif(n = 4^3), dim = c(4, 4, 4))
    edge_probs <- round(edge_probs, 3)
    topic_interaction_patterns <- c(0, 1, 1, 2, 3)
    document_edge_values <- c(0, 0, 1, 1)
    topic <- 1
    tokens_in_document <- 24
    document_sender <- 0
    current_token_topic_assignment <- 1
    current_document_topic_counts <- c(4, 10, 2, 6, 1)

    # first lets try without covariates
    result <- test_internal_functions(
        Test_LSM_Contribution = TRUE,
        edge_probabilities = edge_probs,
        tokens_in_document = tokens_in_document,
        topic = topic - 1,
        current_token_topic_assignment = current_token_topic_assignment,
        current_document_topic_counts = current_document_topic_counts,
        document_edge_values = document_edge_values,
        topic_interaction_patterns = topic_interaction_patterns,
        document_sender = document_sender)

    contribution <- 0
    for (i in 0:(length(document_edge_values) - 1)) {
      if (i != document_sender) {
        sum_term <- test_internal_functions(Test_Sum_Over_T_Edge_Probs = TRUE,
                                            edge_probabilities = edge_probs,
                                            tokens_in_document = tokens_in_document,
                                            current_token_topic_assignment = current_token_topic_assignment,
                                            current_document_topic_counts = current_document_topic_counts,
                                            leave_out_current_token = TRUE,
                                            topic_interaction_patterns = topic_interaction_patterns,
                                            document_sender = document_sender,
                                            document_recipient = i,
                                            leave_out_topic = -1)
        sum_term <- sum_term + 1 / tokens_in_document *
          edge_probs[document_sender + 1, i + 1, topic_interaction_patterns[topic] + 1]
        
        if (document_edge_values[i + 1] == 1) {
          contribution <- contribution + log(sum_term)
        } else {
          contribution <- contribution + log(2 - sum_term)
        }
      }
    }

    expect_equal(contribution, result)
})
