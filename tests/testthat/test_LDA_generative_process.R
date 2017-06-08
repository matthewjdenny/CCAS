test_that("That both LDA generative processes work", {
    skip_on_cran()

    # we need a very simple model with a small number of actors, topics, odcuments etc.
    # create an example distribution
    set.seed(12345)

    resample_token_word_types = FALSE
    num_documents = 5
    words_per_doc = 4
    num_topics = 3
    num_terms = 5
    total_word_count = num_documents * words_per_doc
    # five topics
    alpha = 1
    temp <- runif(num_topics)
    alpha_m = alpha * temp/sum(temp)

    # ten unique words
    beta = 2
    temp2 <- runif(num_terms)
    beta_n = beta * temp2/sum(temp2)

    token_topic_assignments <- vector(mode = "list",
                                      length = num_documents)
    token_word_types <- vector(mode = "list",
                               length = num_documents)

    for (i in 1:num_documents) {
        # initialize to -1 to make sure we throw an error if the C++ generative
        # process does not take care of it.
        token_topic_assignments[[i]] <- rep(-1, words_per_doc)
        # if we are not resampling then they need to be fixed from the start.
        token_word_types[[i]] <- floor(runif(n = words_per_doc,
                                             min = 0,
                                             max = num_terms - .000001))

    }

    traditional <- test_internal_functions(
        Test_Sample_Token_Topics_From_Generative_Process = TRUE,
        alpha_m = alpha_m,
        beta_n = beta_n,
        num_documents = num_documents,
        words_per_doc = words_per_doc,
        token_topic_assignments = token_topic_assignments,
        token_word_types = token_word_types,
        resample_word_types = resample_token_word_types,
        use_collapsed_topic_sampling = FALSE)

    collapsed <- test_internal_functions(
        Test_Sample_Token_Topics_From_Generative_Process = TRUE,
        alpha_m = alpha_m,
        beta_n = beta_n,
        num_documents = num_documents,
        words_per_doc = words_per_doc,
        token_topic_assignments = token_topic_assignments,
        token_word_types = token_word_types,
        resample_word_types = resample_token_word_types,
        use_collapsed_topic_sampling = TRUE)

    traditional2 <- test_internal_functions(
        Test_Sample_Token_Topics_From_Generative_Process = TRUE,
        alpha_m = alpha_m,
        beta_n = beta_n,
        num_documents = num_documents,
        words_per_doc = words_per_doc,
        token_topic_assignments = token_topic_assignments,
        token_word_types = token_word_types,
        resample_word_types = TRUE,
        use_collapsed_topic_sampling = FALSE)

    collapsed2 <- test_internal_functions(
        Test_Sample_Token_Topics_From_Generative_Process = TRUE,
        alpha_m = alpha_m,
        beta_n = beta_n,
        num_documents = num_documents,
        words_per_doc = words_per_doc,
        token_topic_assignments = token_topic_assignments,
        token_word_types = token_word_types,
        resample_word_types = TRUE,
        use_collapsed_topic_sampling = TRUE)

})
