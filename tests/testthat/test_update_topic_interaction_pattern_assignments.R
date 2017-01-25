test_that("That Update_Topic_Interaction_Pattern_Assignments works", {
    skip_on_cran()

    # create an example distribution
    seed = 12345
    set.seed(seed)

    num_documents = 20
    words_per_doc = 10
    num_topics = 5
    num_terms = 10
    num_actors = 4
    num_ip = 4
    numld = 2
    num_covar = 4

    intercepts <- rnorm(num_ip,mean = 0, sd = 2)
    latent_pos <- array(data = rnorm(n = num_ip*num_actors*numld, mean = 0, sd = 2),
                        dim = c(num_ip,num_actors,numld))
    coefficients <- matrix(rnorm(n = num_covar*num_ip, mean = 0, sd = 2),
                           nrow = num_covar, ncol = num_ip)

    document_edge_matrix <- matrix(round(runif(num_actors * num_documents)),
                                   nrow = num_documents,
                                   ncol = num_actors)

    author_indexes <- floor(runif(n = num_documents, min = 0, max = num_actors - 0.00001))

    covars <- array(data = rnorm(n = num_actors*num_actors*num_covar,0,2),
                    dim = c(num_actors, num_actors, num_covar))

    # create token topic assignemnts and token word types list
    init_token_topic_assignments <- vector(mode = "list",
                                           length = num_documents)
    token_word_types <- vector(mode = "list",
                               length = num_documents)
    document_topic_counts <- matrix(0,nrow = num_documents,
                                    ncol = num_topics)
    word_type_topic_counts <- matrix(0,
                                     nrow = num_terms,
                                     ncol = num_topics)
    for (i in 1:num_documents) {
        init_token_topic_assignments[[i]] <- floor(runif(n = words_per_doc,
                                                         min = 0,
                                                         max = num_topics - .000001))
        token_word_types[[i]] <- floor(runif(n = words_per_doc,
                                             min = 0,
                                             max = num_terms - .000001))
        for(j in 0:(num_topics-1)){
            document_topic_counts[i,j] <- length(which(init_token_topic_assignments[[i]]== j))
        }
        for(l in 1:words_per_doc){
            for(j in 1:num_terms){
                for(k in 1:num_topics){
                    if(token_word_types[[i]][l] == (j - 1) & init_token_topic_assignments[[i]][l] == (k -1)){
                        word_type_topic_counts[j,k] <- word_type_topic_counts[j,k] + 1
                    }
                }
            }
        }
    }

    topic_token_counts <- colSums(word_type_topic_counts)

    edge_probs <- array(data = runif(n = 64), dim = c(4,4,4))

    topic_interaction_patterns = c(0,1,1,2,3)
    random_numbers <- runif(num_topics)

    # first lets try without covariates
    result <- test_internal_functions(
        Test_Update_Topic_Interaction_Pattern_Assignments = TRUE,
        author_indexes = author_indexes,
        document_edge_matrix = document_edge_matrix,
        document_topic_counts = document_topic_counts,
        topic_interaction_patterns = topic_interaction_patterns,
        intercepts = intercepts,
        coefficients = coefficients,
        latent_positions = latent_pos,
        covariates = covars,
        using_coefficients = TRUE,
        random_numbers = random_numbers,
        edge_probabilities = edge_probs)

    number_of_documents <- nrow(document_edge_matrix)
    number_of_actors <- ncol(document_edge_matrix)
    number_of_interaction_patterns <- length(intercepts)
    number_of_topics <- length(topic_interaction_patterns)
    number_counter <- 1

    for (t in 1:number_of_topics)  {
      interaction_pattern_assignment_log_probs <- rep(0, number_of_interaction_patterns)
      held_out_sum_over_t_terms <- matrix(0, number_of_documents, number_of_actors)
      for (i in 1:number_of_documents) {
        current_document_topic_counts <- document_topic_counts[i, ]
        document_edge_values <- document_edge_matrix[i, ]
        tokens_in_document <- sum(current_document_topic_counts)
        document_sender <- author_indexes[i]

        for (j in 1:number_of_actors) {
          if (document_sender != (j - 1)) {
            temp <- test_internal_functions(
              Test_Sum_Over_T_Edge_Probs = TRUE,
              edge_probabilities = edge_probs,
              tokens_in_document = tokens_in_document,
              current_token_topic_assignment = -1,
              current_document_topic_counts = current_document_topic_counts,
              leave_out_current_token = FALSE,
              topic_interaction_patterns = topic_interaction_patterns,
              document_sender = document_sender,
              document_recipient = j - 1,
              leave_out_topic = t - 1
            )
            if (document_edge_values[j] == 1)
              held_out_sum_over_t_terms[i, j] <- temp
            else
              held_out_sum_over_t_terms[i, j] <- 1 - temp
          }
        }
      }
      for (c in 1:number_of_interaction_patterns) {
        log_prob <- 0
        for (i in 1:number_of_documents) {
          current_document_topic_counts <- document_topic_counts[i, ]
          document_edge_values <- document_edge_matrix[i, ]
          tokens_in_document <- sum(current_document_topic_counts)
          document_sender <- author_indexes[i]

          for (j in 1:number_of_actors) {
            if (document_sender != (j - 1)) {
              ttc <- current_document_topic_counts[t] / tokens_in_document
              if (document_edge_values[j] == 1) {
                log_prob <- log_prob +
                  log(held_out_sum_over_t_terms[i, j] +
                        ttc * edge_probs[document_sender + 1, j, c])
              } else {
                log_prob <- log_prob +
                  log(held_out_sum_over_t_terms[i, j] +
                        1 - (ttc * edge_probs[document_sender + 1, j, c]))
              }
            }
          }
        }
        interaction_pattern_assignment_log_probs[c] <- log_prob
      }
      random_number <- random_numbers[number_counter]
      number_counter <- number_counter + 1
      new_assignment <- test_internal_functions(
        Test_Log_Space_Multinomial_Sampler = TRUE,
        distribution = interaction_pattern_assignment_log_probs,
        u = random_number, seed = seed
      )
      # see lines 2200-2219: it incrementes by one in c++
      topic_interaction_patterns[t] <- new_assignment - 1
    }

  expect_that(topic_interaction_patterns, equals(as.vector(result)))
})
