test_that("That the sum over t edge probabilities works", {
    skip_on_cran()

    # create an example distribution
    set.seed(12345)

    edge_probs <- array(data = runif(n = 64), dim = c(4,4,4))

    # first lets try without covariates
    result <- test_internal_functions(
        Test_Sum_Over_T_Edge_Probs = TRUE,
        edge_probs = edge_probs,
        tokens_in_document = 24,
        current_token_topic_assignment = 2,
        current_document_topic_counts = c(5,10,2,6,1),
        leave_out_current_token = FALSE,
        topic_interaction_patterns = c(0,1,1,2,3),
        document_sender = 0, # note that these have to be zero indexed
        document_recipient = 2, # note that these have to be zero indexed
        leave_out_topic = -1)


})
