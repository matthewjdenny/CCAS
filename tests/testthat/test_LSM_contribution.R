test_that("That LSM contribution works", {
    skip_on_cran()

    # create an example distribution
    set.seed(12345)

    edge_probs <- array(data = runif(n = 64), dim = c(4,4,4))

    # first lets try without covariates
    result <- test_internal_functions(
        Test_LSM_Contribution = TRUE,
        edge_probabilities = edge_probs,
        tokens_in_document = 24,
        topic = 1,
        current_token_topic_assignment = 2,
        current_document_topic_counts = c(5,10,2,6,1),
        document_edge_values = c(0,0,1,1),
        topic_interaction_patterns = c(0,1,1,2,3),
        document_sender = 0)

    # no errors, will write an analytical test

})
