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
             token_word_type_list = "list",
             token_topic_assignment_list = "list",
             network_covariates_list = "list",
             aggregate_network = "matrix",
             blank_documents = "numeric",
             using_covariates = "logical",
             token_word_type_list_zero_indexed = "list",
             token_topic_assignment_list_zero_indexed = "list"
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
             ComNet_Object = "ComNet",
             MCMC_output = "list",
             parameter_estimates = "data.frame",
             topic_model_results = "list",
             interaction_patterns = "numeric",
             number_of_topics = "numeric",
             latent_space_dimensions = "numeric",
             target_accept_rate = "numeric",
             tollerance = "numeric",
             formula = "formula",
             alpha = "numeric",
             beta = "numeric",
             iterations = "numeric",
             burnin = "numeric",
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
