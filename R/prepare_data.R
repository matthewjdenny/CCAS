#' @title Prepare communication network data for analysis.
#' @description Reads in and stores communication network data in a ComNet
#' object for use in main estimation procedure.
#'
#' @param document_authors A vector recording the index of the sender of each
#' document. This vector must be equal in length to the number of documents
#' in the corpus, and should be one-indexed.
#' @param document_edge_matrix A documents (rows) by actors (columns) matrix
#' recording whether actor j was a recipient of message i. This matrix should
#' be one-indexed.
#' @param document_term_matrix A documents (rows) by unique terms (columns)
#'  matrix recording the count of unique term j in document i.
#' @param covariate_data Defaults to NULL. The user may provide a data.frame
#' containing covariate data for each actor. This data frame must have the same
#' number of rows as actors, and should use descriptive variable names.
#' @param vocabulary Defaults to NULL, in which case the column names of the
#' document_term_matrix will be used as the vocabulary.
#' @return An object of class ComNet that will be used by the ccas() function.
#' @examples
#' \dontrun{
#' # load in example county government email data.
#' data(author_attributes)
#' data(document_edge_matrix)
#' data(document_word_matrix)
#' data(vocabulary)
#'
#' # the first colun of the doc-edge matrix is the author index. Take it out and
#' # then remove it from the doc-edge matrix.
#' document_authors <- document_edge_matrix[,1]
#' document_edge_matrix <- document_edge_matrix[,-1]
#' ComNet <- prepare_data(
#' document_authors = document_authors,
#' document_edge_matrix = document_edge_matrix,
#' document_term_matrix = document_word_matrix,
#' covariate_data = author_attributes,
#' vocabulary = vocabulary)
#' }
#' @export
prepare_data <- function(document_authors,
                         document_edge_matrix,
                         document_term_matrix,
                         covariate_data = NULL,
                         vocabulary = NULL){

    document_edge_matrix <- as.matrix(document_edge_matrix)
    document_term_matrix <- as.matrix(document_term_matrix)
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
    cat("Vocabulary size:", vocabulary_size, "unique terms\n")
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
    token_word_type_list <- vector(mode = "list",
                                   length = num_documents)
    token_topic_assignment_list <- vector(mode = "list",
                                          length = num_documents)
    token_word_type_list_zero_indexed <- vector(mode = "list",
                                                length = num_documents)
    token_topic_assignment_list_zero_indexed <- vector(mode = "list",
                                                       length = num_documents)
    aggregate_network <- matrix(0, nrow = num_actors,ncol = num_actors)
    blank_documents <- rep(0 , num_documents)
    for (i in 1:num_documents) {
        # deal with token_word_type_list, token_topic_assignment_list
        tokens_in_current_doc <- sum(document_term_matrix[i,])

        # if there is atleast one token in the document
        if (tokens_in_current_doc > 0) {
            cur_types <- rep(0, tokens_in_current_doc)
            cur_topics <- rep(0, tokens_in_current_doc)

            # fin out which word types occur a positive nubmer of times in the
            # document
            inds <- which(document_term_matrix[i,] > 0)

            # this counter indicates which entry in the vector we are filling
            counter <- 1

            # loop over each unique word type
            for (j in 1:length(inds)) {
                # loop over the numer of times the term appears in the document
                for (k in 1:document_term_matrix[i,inds[j]]) {
                    # record the index of the token type in the types vector
                    cur_types[counter] <- inds[j]

                    # increment the counter
                    counter <- counter + 1
                }
            }

            # assign the vectors to a list
            token_word_type_list[[i]] <- cur_types
            token_topic_assignment_list[[i]] <- cur_topics
            token_word_type_list_zero_indexed[[i]]  <- cur_types - 1
            # dont need to subtract one because they are all zero to begin with
            token_topic_assignment_list_zero_indexed[[i]] <- cur_topics

        } else {
            # if there were no words in the current document, record this.
            blank_documents[i] <- 1
        }

        # now deal with edge matrix
        cur_edges <- as.numeric(document_edge_matrix[i,])
        if (sum(cur_edges) > 0) {
            aggregate_network[document_authors[i],] <-
                aggregate_network[document_authors[i],] + cur_edges
        } else {
            stop(paste("Document",i,"does not have any recipients, please remove"))
        }

    } # end of loop over documents

    # for now, we do not allow the user to include emails that do not have words
    # in them in the data. Eventually, we will deal with this in the model.
    if (sum(blank_documents) > 0) {
        inds <- which(blank_documents == 1)
        # remove those documents
        document_authors <- document_authors[-inds]
        document_term_matrix <- document_term_matrix[-inds,]
        document_edge_matrix <- document_edge_matrix[-inds,]
        token_word_type_list <- token_word_type_list[-inds]
        token_topic_assignment_list <- token_topic_assignment_list[-inds]
        token_word_type_list_zero_indexed <- token_word_type_list_zero_indexed[-inds]
        token_topic_assignment_list_zero_indexed <- token_topic_assignment_list_zero_indexed[-inds]
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
        blank_documents = blank_documents,
        using_covariates = using_covariates,
        token_word_type_list_zero_indexed = token_word_type_list_zero_indexed,
        token_topic_assignment_list_zero_indexed = token_topic_assignment_list_zero_indexed,
        document_authors_zero_indexed = document_authors - 1)

    # only assign if covariate_data is not NULL in order to aviod error
    if (using_covariates) {
        # make sure that all factor variables are converted to strings
        for (i in 1:ncol(covariate_data)) {
            if (class(covariate_data[,i]) == "factor") {
                covariate_data[,i] <- as.character(covariate_data[,i])
            }
        }

        # now assign to slot in object
        ComNet_Object@covariate_data <- covariate_data
    }

    return(ComNet_Object)
}
