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
#' 'leave_out_topic'. must be provided.
#' @param Test_Prior_Pobability_Of_I_P_Params Defaults to FALSE. If TRUE, then optional
#' arguments intercepts, coefficients, latent_pos, intercept_prior_mean,
#' intercept_prior_variance, coefficient_prior_mean, coefficient_prior_variance,
#' latent_position_prior_mean, latent_position_prior_variance,
#' using_coefficients must be provided.
#' @param Test_Sample_New_I_P_Parameters Defaults to FALSE. If TRUE, then
#' optional  arguments intercepts, coefficients, latent_pos,
#' intercept_proposal_variances, coefficient_proposal_variances,
#' latent_position_proposal_variances, (all of which must be of length = number
#' of interaction patterns, and should all be equal) using_coefficients must be
#' provided.
#' @param Test_LSM_Contribution Defaults to FALSE. If TRUE, then optional
#' arguments edge_probs, tokens_in_document, topic,
#' current_token_topic_assignment, current_document_topic_counts,
#' document_edge_values, topic_interaction_patterns, document_sender must be
#'  provided.
#' @param Test_LDA_Contribution Defaults to FALSE. If TRUE, then optional
#' arguments tokens_in_document, current_token_topic_assignment,
#' current_document_topic_counts, word_type_topic_counts, topic_token_counts,
#' topic, current_word_type, alpha_m, beta_n, beta  must be provided.
#' @param Test_Update_Single_Token_Topic_Assignment Defaults to FALSE. If TRUE,
#' then optional arguments edge_probs, tokens_in_document,
#' current_token_topic_assignment, current_document_topic_counts,
#' word_type_topic_counts, topic_token_counts, current_word_type, alpha_m,
#' beta_n, document_edge_values, topic_interaction_patterns, document_sender,
#' rand_num must be provided.
#' @param Test_Update_All_Token_Topic_Assignments Defaults to FALSE. If TRUE,
#' then optional arguments author_indexes, document_edge_matrix,
#' topic_interaction_patterns, document_topic_counts, word_type_topic_counts,
#' topic_token_counts, token_topic_assignments, token_word_types, intercepts,
#' coefficients, latent_pos, covars, alpha_m, beta_n, random_numbers,
#' using_coefficients must be provided. Make sure that random_numbers has length
#' equal to the total number of tokens in the corpus.
#' @param Test_Update_Interaction_Pattern_Parameters Defaults to FALSE. If TRUE,
#' then optional arguments author_indexes, document_edge_matrix,
#' document_topic_counts, topic_interaction_patterns, intercepts, coefficients,
#' latent_pos, covars, using_coefficients, intercept_prior_mean,
#' intercept_prior_variance, intercept_proposal_variances,
#' coefficient_prior_mean, coefficient_prior_variance,
#' coefficient_proposal_variances, latent_position_prior_mean,
#' latent_position_prior_variance, latent_position_proposal_variances, rand_num,
#' edge_probs must be provided.
#' @param Test_Update_Topic_Interaction_Pattern_Assignments Defaults to FALSE.
#' If TRUE, then optional arguments author_indexes, document_edge_matrix,
#' document_topic_counts, topic_interaction_patterns, intercepts, coefficients,
#' latent_pos, covars, using_coefficients, random_numbers, edge_probs must be
#' provided.
#' @param envir Should not be changed by the user, captures the current
#' environment to facilitate testing.
#' @param ... optional arguments necessary to run each of the internal functions.
#' @return Whatever is returned by the internal function being tested
#' @export
test_internal_functions <- function(
    Test_Log_Space_Multinomial_Sampler = FALSE,
    Test_Edge_Probability = FALSE,
    Test_Sum_Over_T_Edge_Probs = FALSE,
    Test_Prior_Pobability_Of_I_P_Params = FALSE,
    Test_Sample_New_I_P_Parameters = FALSE,
    Test_LSM_Contribution = FALSE,
    Test_LDA_Contribution = FALSE,
    Test_Update_Single_Token_Topic_Assignment = FALSE,
    Test_Update_All_Token_Topic_Assignments = FALSE,
    Test_Update_Interaction_Pattern_Parameters = FALSE,
    Test_Update_Topic_Interaction_Pattern_Assignments = FALSE,
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
        edge_probs <- NULL
        tokens_in_document <- NULL
        current_token_topic_assignment <- NULL
        current_document_topic_counts <- NULL
        leave_out_current_token <- NULL
        topic_interaction_patterns <- NULL
        document_sender <- NULL
        document_recipient <- NULL
        leave_out_topic <- NULL
        intercepts <- NULL
        coefficients <- NULL
        latent_pos <- NULL
        intercept_prior_mean <- NULL
        intercept_prior_variance <- NULL
        intercept_proposal_variances <- NULL
        coefficient_prior_mean <- NULL
        coefficient_prior_variance <- NULL
        coefficient_proposal_variances <- NULL
        latent_position_prior_mean <- NULL
        latent_position_prior_variance <- NULL
        latent_position_proposal_variances <- NULL
        topic <- NULL
        document_edge_values <- NULL
        word_type_topic_counts <- NULL
        topic_token_counts <- NULL
        current_word_type <- NULL
        alpha_m <- NULL
        beta_n <- NULL
        rand_num <- NULL
        author_indexes <- NULL
        document_edge_matrix <- NULL
        document_topic_counts <- NULL
        token_topic_assignments <- NULL
        token_word_types <- NULL
        covars <- NULL
        random_numbers <- NULL
    }

    object <- as.list(substitute(list(...)))[-1L]
    if (length(object) > 0) {
        for (i in 1:length(object)) {
            # try both direct assignment and get()
            if (typeof(object[[i]]) == "symbol") {
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
    if (Test_Log_Space_Multinomial_Sampler) {
        if (is.null(distribution)) {
            stop("you must provide an optional variable 'distribution', which is an unnormalized vetor of log probabilities.")
        }
        if (is.null(seed)) {
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

    if (Test_Prior_Pobability_Of_I_P_Params) {
        return_object <- ppipp(intercepts,
              coefficients,
              latent_pos,
              intercept_prior_mean,
              intercept_prior_variance,
              coefficient_prior_mean,
              coefficient_prior_variance,
              latent_position_prior_mean,
              latent_position_prior_variance,
              using_coefficients)
    }

    if (Test_Sample_New_I_P_Parameters) {
        return_object <- snipp(intercepts,
              coefficients,
              latent_pos,
              intercept_proposal_variances,
              coefficient_proposal_variances,
              latent_position_proposal_variances,
              using_coefficients)
    }

    if (Test_LSM_Contribution) {
        return_object <- lsmc(
            edge_probs,
            tokens_in_document,
            topic,
            current_token_topic_assignment,
            current_document_topic_counts,
            document_edge_values,
            topic_interaction_patterns,
            document_sender)
    }

    if (Test_LDA_Contribution) {
        return_object <- ldac(
            tokens_in_document,
            current_token_topic_assignment,
            current_document_topic_counts,
            word_type_topic_counts,
            topic_token_counts,
            topic,
            current_word_type,
            alpha_m,
            beta_n,
            beta)
    }

    if (Test_Update_Single_Token_Topic_Assignment) {
        return_object <- ustta(edge_probs,
              tokens_in_document,
              current_token_topic_assignment,
              current_document_topic_counts,
              word_type_topic_counts,
              topic_token_counts,
              current_word_type,
              alpha_m,
              beta_n,
              document_edge_values,
              topic_interaction_patterns,
              document_sender,
              rand_num)
    }

    if (Test_Update_All_Token_Topic_Assignments) {
        return_object <-  utta(author_indexes,
            document_edge_matrix,
            topic_interaction_patterns,
            document_topic_counts,
            word_type_topic_counts,
            topic_token_counts,
            token_topic_assignments,
            token_word_types,
            intercepts,
            coefficients,
            latent_pos,
            covars,
            alpha_m,
            beta_n,
            random_numbers,
            using_coefficients)
    }

    if (Test_Update_Interaction_Pattern_Parameters) {
        return_object <- uipp(
            author_indexes,
            document_edge_matrix,
            document_topic_counts,
            topic_interaction_patterns,
            intercepts,
            coefficients,
            latent_pos,
            covars,
            using_coefficients,
            intercept_prior_mean,
            intercept_prior_variance,
            intercept_proposal_variances,
            coefficient_prior_mean,
            coefficient_prior_variance,
            coefficient_proposal_variances,
            latent_position_prior_mean,
            latent_position_prior_variance,
            latent_position_proposal_variances,
            rand_num,
            edge_probs)
    }

    if (Test_Update_Topic_Interaction_Pattern_Assignments) {
        return_object <- utipa(author_indexes,
                               document_edge_matrix,
                               document_topic_counts,
                               topic_interaction_patterns,
                               intercepts,
                               coefficients,
                               latent_pos,
                               covars,
                               using_coefficients,
                               random_numbers,
                               edge_probs)
    }


    # return whatever needs to be returned
    return(return_object)
}
