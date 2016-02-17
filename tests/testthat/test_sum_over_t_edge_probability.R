test_that("That the sum over t edge probabilities works", {
    skip_on_cran()

    # create an example distribution
    set.seed(12345)

    edge_probs <- array(data = runif(n = 64), dim = c(4,4,4))

    # first lets try without holding out a token or topic
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

    # now calculate by hand, 5 topics
    tokens_in_document = 24
    current_document_topic_counts = c(5,10,2,6,1)
    topic_interaction_patterns = c(1,2,2,3,4) #increment in R
    document_sender = 1 #increment in R
    document_recipient = 3 #increment in R
    sum <- 0
    for(i in 1:5){
        sum <- sum + (current_document_topic_counts[i]/tokens_in_document)*
            edge_probs[document_sender,document_recipient,
                       topic_interaction_patterns[i]]
    }

    # check to see if they are equal
    expect_equal(sum,result)

    # now lets try holding out a topic
    result2 <- test_internal_functions(
        Test_Sum_Over_T_Edge_Probs = TRUE,
        edge_probs = edge_probs,
        tokens_in_document = 24,
        current_token_topic_assignment = 2,
        current_document_topic_counts = c(5,10,2,6,1),
        leave_out_current_token = TRUE,
        topic_interaction_patterns = c(0,1,1,2,3),
        document_sender = 0, # note that these have to be zero indexed
        document_recipient = 2, # note that these have to be zero indexed
        leave_out_topic = -1)

    # now calculate by hand, leaving out token in topic 3
    current_document_topic_counts = c(5,10,1,6,1)
    sum2 <- 0
    for(i in 1:5){
        sum2 <- sum2 + (current_document_topic_counts[i]/tokens_in_document)*
            edge_probs[document_sender,document_recipient,
                       topic_interaction_patterns[i]]
    }

    # check to see if they are equal
    expect_equal(sum2,result2)


    # now lets try holding out a topic
    result3 <- test_internal_functions(
        Test_Sum_Over_T_Edge_Probs = TRUE,
        edge_probs = edge_probs,
        tokens_in_document = 24,
        current_token_topic_assignment = 2,
        current_document_topic_counts = c(5,10,2,6,1),
        leave_out_current_token = FALSE,
        topic_interaction_patterns = c(0,1,1,2,3),
        document_sender = 0, # note that these have to be zero indexed
        document_recipient = 2, # note that these have to be zero indexed
        leave_out_topic = 3)

    # now calculate by hand, leaving out token in topic 3
    current_document_topic_counts = c(5,10,2,6,1) # reset
    sum3 <- 0
    for (i in 1:5) {
        if (i != 4) {
            sum3 <- sum3 + (current_document_topic_counts[i]/tokens_in_document)*
                edge_probs[document_sender,document_recipient,
                           topic_interaction_patterns[i]]
        }

    }

    # check to see if they are equal
    expect_equal(sum3,result3)

})
