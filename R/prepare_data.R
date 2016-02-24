#' @title Prepare communication network data for analysis.
#' @description Reads in and stores communication network data in a ComNet
#' object for use in main estimation procedure.
#'
#' @param document_authors A vector recording the index of the sender of each
#' document.
#' @param document_edge_matrix A documents x actors matrix recording whether
#' actor j was a recipient of message i.
#' @param document_term_matrix A documents x unique terms matrix recording the
#' count of unique term j in document i.
#' @param covariate_data Defaults to NULL.
#' @param vocabulary Defaults to NULL.
#' @return An object of class ComNet
prepare_data <- function(document_authors,
                         document_edge_matrix,
                         document_term_matrix,
                         covariate_data = NULL,
                         vocabulary = NULL){

    # initialize boolean indicating whether covariate data was provided
    using_covariates <- FALSE
    if (!is.null(covariate_data)) {
        using_covariates <- TRUE
    }

    # if no vocabulary was provided, get the column names from the
    # document_term_matrix
    if (is.null(vocabulary)) {
        vocabulary <- colnames(document_term_matrix)
    }

    # make sure that the number of documents implied by all input data
    # structures is equal.
    if (nrow(document_edge_matrix) == nrow(document_term_matrix) &
        nrow(document_term_matrix) == length(document_authors)) {
        cat("Number of documents:",length(document_authors),"\n")
    } else {
        stop("nrow(document_edge_matrix), nrow(document_term_matrix), length(document_authors) must all be equal, you provided data objects implying a different number of documents.")
    }

    # calculate some corpus level statistics
    vocabulary_size <- length(vocabulary)
    cat("Vocabulary size:", vocabulary_size, "unique terms...\n")
    num_actors <- ncol(document_edge_matrix)
    cat("Number of actors:", num_actors, "\n")
    num_documents <- nrow(document_edge_matrix)
    num_edges <- sum(document_edge_matrix)
    cat("Total number of message recipients:", num_edges, "\n")

    ComNet_Object <- new("ComNet",
        document_authors = document_authors,
        document_term_matrix = document_term_matrix,
        document_edge_matrix = document_edge_matrix,
        vocabulary = vocabulary,
        num_tokens = num_tokens,
        num_actors = num_actors,
        num_documents = num_documents,
        num_edges = num_edges,
        vocabulary_size = vocabulary_size,
        token_word_type_list = token_word_type_list,
        token_topic_assignment_list = token_topic_assignment_list,
        network_covariates_list = network_covariates_list,
        aggregate_network = aggregate_network)

    # only assign if covariate_data is not NULL in order to aviod error
    if (using_covariates) {
        ComNet_Object@covariate_data <- covariate_data
    }

    return(ComNet_Object)
}
