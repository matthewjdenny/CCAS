#' @title Test C++ functions
#' @description Allows the user to interface with the individual C++ functions.
#'
#' @param Test_Log_Space_Multinomial_Sampler Defualts to FALSE. If TRUE, then
#' optional arguments 'distribution' and 'seed' must be provided.
#' @param Test_Edge_Probability Defaults to FALSE. If TRUE, then optional arguments
#' 'intercepts', 'coefficients', 'latent_pos', 'sender', 'recipient',
#' 'current_covariates', 'interaction_pattern', and 'using_coefficients' must be
#' provided.
#' @param Test_Sum_Over_T_Edge_Probs Defaults to FALSE. If TRUE, then optional
#' arguments 'edge_probs', 'tokens_in_document', 'current_token_topic_assignment',
#' 'current_document_topic_counts', 'leave_out_current_token',
#' 'topic_interaction_patterns', 'document_sender','document_recipient', and
#' 'leave_out_topic'.
#' @param envir Should not be changed by the user, captures the current environment
#' to facilitate testing.
#' @param ... optional arguments necessary to run each of the internal functions.
#' @return Whatever is returned by the internal function being tested
#' @export
test_internal_functions <- function(Test_Log_Space_Multinomial_Sampler = FALSE,
                                    Test_Edge_Probability = FALSE,
                                    Test_Sum_Over_T_Edge_Probs = FALSE,
                                    envir = environment(),
                                    ...){

    # set all optional variables potentially used in the function to NULL inside
    # a conditional that is never satisfied so that we avoid R CMD check error.
    dont_preassign <- TRUE
    if(!dont_preassign){
        return_object <- NULL
        seed <- NULL
        distribution <- NULL
        intercepts <- NULL
        coefficients <- NULL
        latent_pos <- NULL
        sender <- NULL
        recipient <- NULL
        current_covariates <- NULL
        interaction_pattern <- NULL
        using_coefficients <- NULL
    }

    object <- as.list(substitute(list(...)))[-1L]
    if (length(object) > 0) {
        for (i in 1:length(object)){
            # try both direct assignment and get()
            if(typeof(object[[i]]) == "symbol"){
                # have to do this double get trick to make it work in all contexts.
                temp <- dynGet(as.character(object[[i]]),
                                       ifnotfound = get(as.character(object[[i]]),
                                                        envir = envir))
                assign(names(object)[i],temp)
            }else{
                assign(names(object)[i],object[[i]])
            }
        }
    }

    # test the lsms function
    if(Test_Log_Space_Multinomial_Sampler){
        if(is.null(distribution)){
            stop("you must provide an optional variable 'distribution', which is an unnormalized vetor of log probabilities.")
        }
        if(is.null(seed)){
            stop("you must provide an optional variable 'seed', which is an integer.")
        }
        return_object <- lsms(distribution,
                              seed)
    }

    # test ep function
    if (Test_Edge_Probability) {
        return_object <- ep(intercepts,
                         coefficients,
                         latent_pos,
                         sender,
                         recipient,
                         current_covariates,
                         interaction_pattern,
                         using_coefficients)
    }

    if (Test_Sum_Over_T_Edge_Probs) {
        return_object <- sotep(edge_probs,
              tokens_in_document,
              current_token_topic_assignment,
              current_document_topic_counts,
              leave_out_current_token,
              topic_interaction_patterns,
              document_sender,
              document_recipient,
              leave_out_topic)
    }


    # return whatever needs to be returned
    return(return_object)
}
