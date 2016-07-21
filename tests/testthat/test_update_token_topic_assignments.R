test_that("That Update_Token_Topic_Assignments works", {
    skip_on_cran()

    # create an example distribution
    seed <- 12345
    set.seed(seed)
    num_documents <- 20
    words_per_doc <- 10
    num_topics <- 5
    num_terms <- 10
    num_actors <- 4
    num_ip <- 4
    numld <- 2
    num_covar <- 4
    total_word_count <- num_documents * words_per_doc

    # five topics
    alpha <- 1
    random_numbers <- runif(total_word_count)
    temp <- runif(num_topics)
    alpha_m <- alpha * temp / sum(temp)

    # ten unique words
    beta <- 2
    temp2 <- runif(num_terms)
    beta_n <- beta * temp2 / sum(temp2)

    intercepts <- rnorm(num_ip,mean = 0, sd = 2)
    latent_pos <- array(data = rnorm(n = num_ip*num_actors*numld, mean = 0, sd = 2),
                        dim = c(num_ip,num_actors,numld))
    coefficients <- matrix(rnorm(n = num_covar*num_ip, mean = 0, sd = 2),
                           nrow = num_covar, ncol = num_ip)

    document_edge_matrix <- matrix(round(runif(num_actors * num_documents)), nrow = num_documents,
                                   ncol = num_actors)

    author_indexes <- floor(runif(n = num_documents, min = 0, max = num_actors - 0.00001))

    covars <- array(data = rnorm(n = num_actors*num_actors * num_covar, 0, 2),
                    dim = c(num_actors, num_actors, num_covar))

    # create token topic assignemnts and token word types list
    init_token_topic_assignments <- vector(mode = "list", length = num_documents)
    save_token_topic_assignments <- vector(mode = "list", length = num_documents)
    token_word_types <- vector(mode = "list", length = num_documents)
    document_topic_counts <- matrix(0, nrow = num_documents, ncol = num_topics)
    word_type_topic_counts <- matrix(0, nrow = num_terms, ncol = num_topics)

    for (i in 1:num_documents) {
        temp <- floor(runif(n = words_per_doc, min = 0, max = num_topics - .000001))
        init_token_topic_assignments[[i]] <- temp
        save_token_topic_assignments[[i]] <- temp
        token_word_types[[i]] <- floor(runif(n = words_per_doc, min = 0, max = num_terms - .000001))
        for(j in 0:(num_topics - 1)) {
          document_topic_counts[i, j] <- length(which(init_token_topic_assignments[[i]] == j))
        }
        for (l in 1:words_per_doc) {
            for (j in 1:num_terms) {
                for (k in 1:num_topics) {
                    if (token_word_types[[i]][l] == (j - 1) & init_token_topic_assignments[[i]][l] == (k - 1)) {
                        word_type_topic_counts[j, k] <- word_type_topic_counts[j, k] + 1
                    }
                }
            }
        }
    }

    number_of_documents <- nrow(document_edge_matrix)
    number_of_actors <- ncol(document_edge_matrix)
    number_of_interaction_patterns <- length(intercepts)
    number_of_topics <- length(alpha_m)
    topic_token_counts <- colSums(word_type_topic_counts)
    token_topic_assignments <- init_token_topic_assignments
    topic_interaction_patterns <- c(0, 1, 1, 2, 3)

    # do.call(cbind, token_topic_assignments)

    # first lets try without covariates
    result <- test_internal_functions(
        Test_Update_All_Token_Topic_Assignments = TRUE,
        author_indexes = author_indexes,
        document_edge_matrix = document_edge_matrix,
        topic_interaction_patterns = topic_interaction_patterns,
        document_topic_counts = document_topic_counts,
        word_type_topic_counts = word_type_topic_counts,
        topic_token_counts = topic_token_counts,
        token_topic_assignments = token_topic_assignments,
        token_word_types = token_word_types,
        intercepts = intercepts,
        coefficients = coefficients,
        latent_positions = latent_pos,
        covariates = covars,
        alpha_m = alpha_m,
        beta_n = beta_n,
        random_numbers = random_numbers,
        using_coefficients = FALSE,
        parallel = FALSE)

    # check to see if anything changed
    # THe strange thing about Rcpp shallow data structures is that they
    # change the original along with the underlying data, so these two things
    # are now equal. But this is not right! now we have the initial and updated
    # token topic assignments equal to eachother. They are even equal to
    # init_token_topic_assignments if you want to check, which is even more
    # interesting since we assigned from it on line 73.
    expect_that(do.call(cbind, result[[4]]), equals(do.call(cbind, token_topic_assignments)))

    # but now noticing that the ending token topic assignments are equal to the
    # starting ones, we need to reset the token_topic_assignments to
    # save_token_topic_assignments (which broke referencing -- see line 52).
    # This puts everything right again. You can check by explicitly printing
    # the various data structures before and after calling the function on line
    # 79.
    # do.call(cbind, save_token_topic_assignments)

    token_topic_assignments <- save_token_topic_assignments

    edge_probabilities <- array(0, c(ncol(document_edge_matrix), ncol(document_edge_matrix), length(intercepts)))
    for (i in 1:ncol(document_edge_matrix)) {
      for (j in 1:ncol(document_edge_matrix)) {
        if (i != j) {
          for (k in 1:length(intercepts))
            edge_probabilities[i, j, k] <- test_internal_functions(
              Test_Edge_Probability = TRUE,
              intercepts = intercepts,
              coefficients = coefficients,
              latent_positions = latent_pos,
              document_sender = i - 1,
              document_recipient = j - 1,
              current_covariates = covars[i, j, ],
              interaction_pattern = k - 1,
              using_coefficients = FALSE
            )
        }
      }
    }
    expect_that(result[[5]], equals(edge_probabilities))

    rand_num_counter <- 1
    for (d in 1:number_of_documents) {
      current_token_topic_assignments <- token_topic_assignments[[d]]
      current_token_word_types <- token_word_types[[d]]
      tokens_in_document <- length(current_token_topic_assignments)
      current_document_topic_counts <- document_topic_counts[d, ]
      document_edge_values <- document_edge_matrix[d, ]
      document_sender <- author_indexes[d]

      for (i in 1:tokens_in_document) {
        current_token_topic_assignment <- current_token_topic_assignments[i]
        current_word_type <- current_token_word_types[i]
        rand_num <- random_numbers[rand_num_counter]
        rand_num_counter <- rand_num_counter + 1

        new_topic_assignment_vec <- test_internal_functions(
          Test_Update_Single_Token_Topic_Assignment = TRUE,
          edge_probabilities = edge_probabilities,
          tokens_in_document = tokens_in_document,
          current_token_topic_assignment = current_token_topic_assignment,
          current_document_topic_counts = current_document_topic_counts,
          word_type_topic_counts = word_type_topic_counts,
          topic_token_counts = topic_token_counts,
          current_word_type = current_word_type,
          alpha_m = alpha_m,
          beta_n = beta_n,
          document_edge_values = document_edge_values,
          topic_interaction_patterns = topic_interaction_patterns,
          document_sender = document_sender,
          rand_num = rand_num,
          parallel = FALSE,
          use_cached_token_topic_distribution = FALSE,
          cached_token_topic_distribution = rep(0, number_of_topics)
        )

        new_topic_assignment <- new_topic_assignment_vec[1]
        if (new_topic_assignment != current_token_topic_assignment) {
          document_topic_counts[d, current_token_topic_assignment + 1] <-
            document_topic_counts[d, current_token_topic_assignment + 1] - 1
          document_topic_counts[d, new_topic_assignment + 1] <-
            document_topic_counts[d, new_topic_assignment + 1] + 1

          current_token_topic_assignments[i] <- new_topic_assignment

          topic_token_counts[current_token_topic_assignment + 1] <-
            topic_token_counts[current_token_topic_assignment + 1] - 1
          topic_token_counts[new_topic_assignment + 1] <-
            topic_token_counts[new_topic_assignment + 1] + 1

          word_type_topic_counts[current_word_type + 1, current_token_topic_assignment + 1] <-
            word_type_topic_counts[current_word_type + 1, current_token_topic_assignment + 1] - 1
          word_type_topic_counts[current_word_type + 1, new_topic_assignment + 1] <-
            word_type_topic_counts[current_word_type + 1, new_topic_assignment + 1] + 1
        }
      }
      token_topic_assignments[[d]] <- current_token_topic_assignments
    }

    # substracting one gives values of negative 1 in
    # do.call(cbind, token_topic_assignments) - 1
    # expect_that(do.call(cbind, result[[4]]),
    # equals(do.call(cbind, token_topic_assignments) - 1))

    # now everything works out!
    expect_that(do.call(cbind, result[[4]]), equals(do.call(cbind, token_topic_assignments)))

})
