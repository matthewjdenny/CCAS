test_that("That Update_Interaction_Pattern_Parameters works", {
    skip_on_cran()

    # create an example distribution
    set.seed(12345)
    rand_num <- runif(1)

    num_documents <- 20
    words_per_doc <- 10
    num_topics <- 5
    num_terms <- 10
    num_actors <- 4
    num_ip <- 4
    numld <- 2
    num_covar <- 4

    intercepts <- rnorm(num_ip, mean = 0, sd = 2)
    latent_pos <- array(data = rnorm(n = num_ip * num_actors * numld, mean = 0, sd = 2),
                        dim = c(num_ip, num_actors, numld))
    coefficients <- matrix(rnorm(n = num_covar * num_ip, mean = 0, sd = 2),
                           nrow = num_covar, ncol = num_ip)
    document_edge_matrix <- matrix(round(runif(num_actors * num_documents)),
                                   nrow = num_documents,
                                   ncol = num_actors)
    author_indexes <- floor(runif(n = num_documents, min = 0, max = num_actors - 0.00001))
    covars <- array(data = rnorm(n = num_actors * num_actors * num_covar, 0, 2),
                    dim = c(num_actors, num_actors, num_covar))

    # create token topic assignemnts and token word types list
    init_token_topic_assignments <- vector(mode = "list", length = num_documents)
    token_word_types <- vector(mode = "list", length = num_documents)
    document_topic_counts <- matrix(0, nrow = num_documents, ncol = num_topics)
    word_type_topic_counts <- matrix(0, nrow = num_terms, ncol = num_topics)
    
    for (i in 1:num_documents) {
      init_token_topic_assignments[[i]] <- floor(runif(n = words_per_doc,
                                                       min = 0,
                                                       max = num_topics - .000001))
      token_word_types[[i]] <- floor(runif(n = words_per_doc,
                                           min = 0,
                                           max = num_terms - .000001))
      for (j in 0:(num_topics - 1)) {
        document_topic_counts[i, j] <- length(which(init_token_topic_assignments[[i]]== j))
      }
      for (l in 1:words_per_doc) {
        for (j in 1:num_terms) {
          for (k in 1:num_topics) {
            if (token_word_types[[i]][l] == (j - 1) & init_token_topic_assignments[[i]][l] == (k - 1)){
              word_type_topic_counts[j, k] <- word_type_topic_counts[j, k] + 1
            }
          }
        }
      }
    }

    topic_token_counts <- colSums(word_type_topic_counts)
    edge_probs <- array(data = runif(n = 64), dim = c(4, 4, 4))
    topic_interaction_patterns = c(0, 1, 1, 2, 3)

    set.seed(12345)
    result <- test_internal_functions(
        Test_Update_Interaction_Pattern_Parameters = TRUE,
        author_indexes = author_indexes,
        document_edge_matrix = document_edge_matrix,
        document_topic_counts = document_topic_counts,
        topic_interaction_patterns = topic_interaction_patterns,
        intercepts = intercepts,
        coefficients = coefficients,
        latent_positions = latent_pos,
        covariates = covars,
        using_coefficients = TRUE,
        intercept_prior_mean = 0,
        intercept_prior_standard_deviation = 4,
        intercept_proposal_standard_deviations = c(0.5, 0.5, 0.5, 0.5),
        coefficient_prior_mean = 0,
        coefficient_prior_standard_deviation = 4,
        coefficient_proposal_standard_deviations = c(0.5, 0.5, 0.5, 0.5),
        latent_position_prior_mean = 0,
        latent_position_prior_standard_deviation = 4,
        latent_position_proposal_standard_deviations = c(0.5, 0.5, 0.5, 0.5),
        rand_num = rand_num,
        edge_probabilities = edge_probs)

    number_of_documents <- nrow(document_edge_matrix)
    number_of_actors <- ncol(document_edge_matrix)
    number_of_interaction_patterns <- length(intercepts)

    set.seed(12345)
    proposed_parameters <- test_internal_functions(
      Test_Sample_New_I_P_Parameters = TRUE,
      intercepts = intercepts,
      coefficients = coefficients,
      latent_positions = latent_pos,
      intercept_proposal_standard_deviations = c(.5, .5, .5, .5),
      coefficient_proposal_standard_deviations = c(.5, .5, .5, .5),
      latent_position_proposal_standard_deviations = c(.5, .5, .5, .5),
      using_coefficients = TRUE
    )

    proposed_intercepts <- proposed_parameters[[1]]
    proposed_coefficients <- proposed_parameters[[2]]
    proposed_latent_positions <- proposed_parameters[[3]]

    log_current_prior <- test_internal_functions(
      Test_Prior_Pobability_Of_I_P_Params = TRUE,
      intercepts = intercepts,
      coefficients = coefficients,
      latent_positions = latent_pos,
      intercept_prior_mean = 0,
      intercept_prior_standard_deviation = 4,
      coefficient_prior_mean = 0,
      coefficient_prior_standard_deviation = 4,
      latent_position_prior_mean = 0,
      latent_position_prior_standard_deviation = 4,
      using_coefficients = TRUE
    )

    log_proposed_prior <- test_internal_functions(
      Test_Prior_Pobability_Of_I_P_Params = TRUE,
      intercepts = proposed_intercepts,
      coefficients = proposed_coefficients,
      latent_positions = proposed_latent_positions,
      intercept_prior_mean = 0,
      intercept_prior_standard_deviation = 4,
      coefficient_prior_mean = 0,
      coefficient_prior_standard_deviation = 4,
      latent_position_prior_mean = 0,
      latent_position_prior_standard_deviation = 4,
      using_coefficients = TRUE
    )

    log_current_probability <- 0
    log_proposed_probability <- 0
    proposed_edge_probabilities <- array(0, c(number_of_actors, number_of_actors, number_of_interaction_patterns))

    for (i in 1:number_of_actors) {
      for (j in 1:number_of_actors) {
        current_covariates <- covars[i, j, ]
        if (i != j) {
          for (k in 1:number_of_interaction_patterns) {
            proposed_edge_probabilities[i, j, k] <- test_internal_functions(
              Test_Edge_Probability = TRUE,
              intercepts = proposed_intercepts,
              coefficients = proposed_coefficients,
              latent_positions = proposed_latent_positions,
              document_sender = i - 1,
              document_recipient = j - 1,
              current_covariates = current_covariates,
              interaction_pattern = k - 1,
              using_coefficients = TRUE
            )
          }
        }
      }
    }

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
            leave_out_topic = -1
          )
          
          temp2 <- test_internal_functions(
            Test_Sum_Over_T_Edge_Probs = TRUE,
            edge_probabilities = proposed_edge_probabilities,
            tokens_in_document = tokens_in_document,
            current_token_topic_assignment = -1,
            current_document_topic_counts = current_document_topic_counts,
            leave_out_current_token = FALSE,
            topic_interaction_patterns = topic_interaction_patterns,
            document_sender = document_sender,
            document_recipient = j - 1,
            leave_out_topic = -1
          )

          if (document_edge_values[j] == 1) {
            log_current_probability <- log_current_probability + log(temp)
            log_proposed_probability <- log_proposed_probability + log(temp2)
          } else {
            log_current_probability <- log_current_probability + log(1 - temp)
            log_proposed_probability <- log_proposed_probability + log(1 - temp2)
          }
        }
      }
    }

    accept_log_prob <- log_proposed_probability + log_proposed_prior -
      log_current_probability - log_current_prior
    log_random_number <- log(rand_num)

    if (accept_log_prob > log_random_number) {
      expect_that(result[[5]], equals(1))
      expect_that(result[[1]], equals(as.matrix(proposed_intercepts)))
      expect_that(result[[2]], equals(proposed_coefficients))
      expect_that(result[[3]], equals(proposed_latent_positions))
      expect_that(result[[4]], equals(proposed_edge_probabilities))
      expect_that(result[[6]], equals(log_current_probability))
      expect_that(result[[7]], equals(intercepts))
      expect_that(result[[8]], equals(coefficients))
      expect_that(result[[9]], equals(latent_positions))
      expect_that(result[[10]], equals(edge_probs))
      expect_that(result[[11]], equals(log_proposed_probability))
    } else {
      expect_that(result[[5]], equals(0))
      expect_that(result[[1]], equals(as.matrix(intercepts)))
      expect_that(result[[2]], equals(coefficients))
      expect_that(result[[3]], equals(latent_pos))
      expect_that(result[[4]], equals(edge_probs))
      expect_that(result[[6]], equals(log_current_probability))
      expect_that(result[[7]], equals(proposed_intercepts))
      expect_that(result[[8]], equals(proposed_coefficients))
      expect_that(result[[9]], equals(proposed_latent_positions))
      expect_that(result[[10]], equals(proposed_edge_probabilities))
      expect_that(result[[11]], equals(log_proposed_probability))
    }
})
