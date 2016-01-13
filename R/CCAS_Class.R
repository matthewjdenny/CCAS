# an S4 class for ccas objects
setClass(Class = "ccas",
         representation = representation(
             document_term_matrix = "matrix",
             document_edge_matrix = "matrix",
             covariate_data = "data.frame",
             vocabulary = "data.frame",
             num_nodes = "numeric",
             MCMC_output = "list",
             parameter_estimates = "data.frame",
             topic_model_results = "list",
             number_of_interaction_patterns = "numeric",
             number_of_topics = "numeric",
             latent_space_dimensions = "numeric",
             formula = "formula",
             alpha = "numeric",
             beta = "numeric",
             iterations = "numeric",
             burnin = "numeric",
             LSM_prior_variance = "numeric"
         ),
         validity = function(object) {
#              if (!"matrix" %in% class(object@network) & is.null(object@network)
#                  == FALSE) {
#                  stop("'network' must be a 'matrix' object or 'NULL'.")
#              }
             return(TRUE)
         }
)

# define coef for pretty output of ccas object
# setMethod(f = "coef", signature = "ccas", definition = function(object, ...) {
#     return(list(Theta  = object@theta.coef, Lambda = object@lambda.coef))
# }
# )

# define 'show' to get a pretty output of the ccas object
setMethod(f = "show", signature = "ccas", definition = function(object){
    message("Number of Nodes:")
    print(object@num_nodes)
    message("Vocabulary Size:")
    print(length(object@vocabulary))
    message("Number of Messages:")
    print(nrow(object@document_edge_matrix))
    message("Parameter Estimates:")
    print(object@parameter_estimates)
}
)
