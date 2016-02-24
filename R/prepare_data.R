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

    # make sure that the author indexes match up to the edge matrix
    if (max(document_authors) > ncol(document_edge_matrix)) {
        stop(paste("The values in document_authors must correspond to the columns in document_edge_matrix. You provided a document_author index that was larger than the number of columns in document_edge_matrix. Please respecify..."))
    }

    # calculate some corpus level statistics
    vocabulary_size <- length(vocabulary)
    cat("Vocabulary size:", vocabulary_size, "unique terms...\n")
    num_actors <- ncol(document_edge_matrix)
    cat("Number of actors:", num_actors, "\n")
    num_documents <- nrow(document_edge_matrix)
    num_edges <- sum(document_edge_matrix)
    cat("Total number of message recipients:", num_edges, "\n")
    num_tokens <- sum(document_term_matrix)
    cat("Total number of tokens in corpus:", num_tokens, "\n")

    # generate token_word_type_list and a blank token_topic_assignment_list that
    # will be filled in later. Also generate an aggregate_network by agregating
    # over the document_edge_matrix

    # first allocated blank data objects
    token_word_type_list <- vector(mode = "list", length = num_documents)
    token_topic_assignment_list <- vector(mode = "list", length = num_documents)
    aggregate_network <- matrix(0, nrow = num_documents, ncol = num_documents)
    blank_documents <- rep(0 , num_documents)
    for (i in 1:num_documents) {
        # deal with token_word_type_list, token_topic_assignment_list
        tokens_in_current_doc <- sum(document_term_matrix[i,])
        cur_types <- rep(0, tokens_in_current_doc)
        cur_topics <- rep(0, tokens_in_current_doc)

        # fin out which word types occur a positive nubmer of times in the
        # document
        inds <- which(document_term_matrix[i,] > 0)

        if (length(inds) > 0) {

        } else {
            # if there were no words in the current document, record this.
            blank_documents[i] <- 1
        }
    }

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
        aggregate_network = aggregate_network,
        blank_documents = blank_documents)

    # only assign if covariate_data is not NULL in order to aviod error
    if (using_covariates) {
        ComNet_Object@covariate_data <- covariate_data
    }

    # network_covariates_list will be left NULL for now but will eventually
    # store any network covariates passed in to the ccas function.

    return(ComNet_Object)
}
