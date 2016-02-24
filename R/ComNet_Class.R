# an S4 class for ccas objects
setClass(Class = "ComNet",
         representation = representation(
             document_authors = "numeric",
             document_term_matrix = "matrix",
             document_edge_matrix = "matrix",
             covariate_data = "data.frame",
             vocabulary = "data.frame",
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
             using_covariates = "logical"
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
