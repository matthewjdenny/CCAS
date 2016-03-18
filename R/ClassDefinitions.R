# an S4 class for ccas objects
setClass(Class = "ComNet",
         representation = representation(
             document_authors = "numeric",
             document_term_matrix = "matrix",
             document_edge_matrix = "matrix",
             covariate_data = "data.frame",
             vocabulary = "character",
             num_tokens = "numeric",
             num_actors = "numeric",
             num_documents = "numeric",
             num_edges = "numeric",
             vocabulary_size = "numeric",
             aggregate_network = "matrix",
             blank_documents = "numeric",
             using_covariates = "logical",
             document_authors_zero_indexed = "numeric",
             token_word_type_list_zero_indexed = "list",
             token_topic_assignment_list_zero_indexed = "list",
             token_word_type_list = "list",
             token_topic_assignment_list = "list",
             network_covariates_list = "list"
         ),
         validity = function(object) {
             #              if (!"matrix" %in% class(object@network) & is.null(object@network)
             #                  == FALSE) {
             #                  stop("'network' must be a 'matrix' object or 'NULL'.")
             #              }
             return(TRUE)
         }
)

# an S4 class for ccas objects
setClass(Class = "CCAS",
         representation = representation(
             parameter_estimates = "data.frame",
             interaction_patterns = "numeric",
             number_of_topics = "numeric",
             latent_space_dimensions = "numeric",
             target_accept_rate = "numeric",
             tollerance = "numeric",
             formula = "formula",
             alpha = "numeric",
             beta = "numeric",
             alpha_m = "numeric",
             beta_n = "numeric",
             iterations = "numeric",
             metropolis_hastings_iterations = "numeric",
             final_metropolis_hastings_iterations = "numeric",
             final_metropolis_hastings_burnin = "numeric",
             thin = "numeric",
             LSM_intercept_proposal_variance = "numeric",
             LSM_intercept_prior_variance = "numeric",
             LSM_intercept_prior_mean = "numeric",
             LSM_position_proposal_variance = "numeric",
             LSM_position_prior_variance = "numeric",
             LSM_position_prior_mean = "numeric",
             LSM_coefficient_proposal_variance = "numeric",
             LSM_coefficient_prior_variance = "numeric",
             LSM_coefficient_prior_mean = "numeric",
             number_of_covariates = "numeric",
             iterations_before_t_i_p_updates = "numeric",
             update_t_i_p_every_x_iterations = "numeric",
             perform_adaptive_metropolis = "logical",
             adaptive_metropolis_update_size = "numeric",
             seed = "numeric",
             covariate_array = "array",
             ComNet_Object = "ComNet",
             MCMC_output = "list",
             topic_model_results = "list",
             latent_variables = "list"
         ),
         validity = function(object) {
             #              if (!"matrix" %in% class(object@network) & is.null(object@network)
             #                  == FALSE) {
             #                  stop("'network' must be a 'matrix' object or 'NULL'.")
             #              }
             return(TRUE)
         }
)
