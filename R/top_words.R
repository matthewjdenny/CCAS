#' @title Generate Topic top words
#' @description Generates data.frames containing the top words associated with
#' each topic and with each interaction pattern
#'
#' @param CCAS_Object The object returned by the ccas() main estimation
#' function.
#' @param output_LaTeX_code Logical indicating whether LaTeX code should be
#' output by the functon for incorporating in tables. Defaults to FALSE.
#' @param top_words_to_print Defualts to 10, cna be adjusted by the user.
#' @return A CCAS object containing a list of data.frames with topic and
#' interaction pattern top words in the 'at' model_output$top_words field.
#' @export
top_words <- function(CCAS_Object,
                      output_LaTeX_code = FALSE,
                      top_words_to_print = 10) {

    # allocate some variables we will use below:
    interaction_patterns <- CCAS_Object@interaction_patterns
    topics <- CCAS_Object@number_of_topics
    vocab_size <- CCAS_Object@ComNet_Object@vocabulary_size
    vocab <- CCAS_Object@ComNet_Object@vocabulary
    topic_ips <- CCAS_Object@topic_model_results$topic_interaction_patterns
    topic_ips <- topic_ips[nrow(topic_ips),]
    topic_word_counts <- CCAS_Object@topic_model_results$word_type_topic_counts
    topic_token_counts <- as.numeric(CCAS_Object@topic_model_results$topic_token_counts)

    # allocate a dataframe for top words and another for ordered counts
    topic_top_words <- matrix(NA, nrow = topics, ncol = vocab_size)
    topic_top_word_counts <- matrix(0, nrow = topics, ncol = vocab_size)

    # loop through topics and generate ordered counts
    for (i in 1:topics) {
        cur_counts <- topic_word_counts[,i]
        # if there were any words assigned to the topic
        if (sum(cur_counts) > 0) {
            ordering <- order(cur_counts, decreasing = TRUE)
            keep <- which(cur_counts[ordering] > 0)
            topic_top_words[i,keep] <- vocab[ordering[keep]]
            topic_top_word_counts[i,] <- cur_counts[ordering]
        }
    }

    # get column and row names right
    colnames(topic_top_words) <- colnames(topic_top_word_counts) <-
        paste("Top_Word_",1:vocab_size, sep = "")
    rownames(topic_top_words) <- rownames(topic_top_word_counts) <-
        paste("Topic_",1:topics, sep = "")

    cat(" ################\n",
        "Topic Top Words:\n",
        "################\n\n")

    for (i in 1:topics) {
        # get the top ten words
        temp <- topic_top_words[i,1:top_words_to_print]
        # remove any NA words
        rem <- which(is.na(temp))
        if (length(rem) > 0) {
            temp <- temp[-rem]
        }
        # print them out
        cat("Topic:",i,"--",temp, "\n")
    }

    # now we get interaction pattern top words
    # allocate a dataframe for top words and another for ordered counts
    interaction_pattern_top_words <- matrix(NA,
        nrow = interaction_patterns, ncol = vocab_size)
    interaction_pattern_top_word_counts <- matrix(0,
        nrow = interaction_patterns, ncol = vocab_size)
    iteraction_pattern_word_counts <- rep(0,interaction_patterns)

    # loop through topics and generate ordered counts
    for (i in 1:interaction_patterns) {
        cols <- which(topic_ips == i)
        if (length(cols) > 0) {
            iteraction_pattern_word_counts[i] = sum(topic_token_counts[cols])
            cur_counts <- apply(topic_word_counts[,cols], 1, sum)
            # if there were any words assigned to the topic
            if (sum(cur_counts) > 0) {
                ordering <- order(cur_counts, decreasing = TRUE)
                keep <- which(cur_counts[ordering] > 0)
                interaction_pattern_top_words[i,keep] <- vocab[ordering[keep]]
                interaction_pattern_top_word_counts[i,] <- cur_counts[ordering]
            }
        }
    }

    # get column and row names right
    colnames(interaction_pattern_top_words) <-
        colnames(interaction_pattern_top_word_counts) <-
        paste("Top_Word_",1:vocab_size, sep = "")
    rownames(interaction_pattern_top_words) <-
        rownames(interaction_pattern_top_word_counts) <-
        paste("Interaction_Pattern_",1:interaction_patterns, sep = "")

    cat("\n\n ##############################\n",
        "Interaction Pattern Top Words:\n",
        "##############################\n\n")
    for (i in 1:interaction_patterns) {
        # get the top ten words
        temp <- interaction_pattern_top_words[i,1:top_words_to_print]
        # remove any NA words
        rem <- which(is.na(temp))
        if (length(rem) > 0) {
            temp <- temp[-rem]
        }
        # print them out
        cat("Interaction Pattern:",i,"--",temp, "\n")
    }

    # put verything in a list
    top_list <- list(
        interaction_pattern_top_words = interaction_pattern_top_words,
        interaction_pattern_top_word_counts = interaction_pattern_top_word_counts,
        topic_top_words = topic_top_words,
        topic_top_word_counts = topic_top_word_counts)
    # put it in the CCAS)Object
    CCAS_Object@model_output$top_words <- top_list

    # now generate LaTeX output if requested:
    if (output_LaTeX_code) {
        word_table <- cbind(topic_token_counts,
                            topic_top_words)
        word_table <- as.data.frame(word_table)

        cat("\n\n #############################\n",
            "Topic Top Words LaTeX Output:\n",
            "#############################\n\n")

        SpeedReader::color_word_table(word_table,
                                     covariate_columns  = 1,
                                     word_columns = 2:(top_words_to_print+1),
                                     min_black = 30,
                                     print_first = topics,
                                     all_same = FALSE,
                                     remove_words = "",
                                     second_table = NULL,
                                     max_char_width = 60,
                                     bold_covariates = T)

        word_table <- cbind(iteraction_pattern_word_counts,
                            interaction_pattern_top_words)
        word_table <- as.data.frame(word_table)

        cat("\n\n ##########################################\n",
            "Interaction Pattern Top Words LaTeX Output:\n",
            "###########################################\n\n")

        SpeedReader::color_word_table(word_table,
                                      covariate_columns  = 1,
                                      word_columns = 2:(top_words_to_print+1),
                                      min_black = 30,
                                      print_first = interaction_patterns,
                                      all_same = FALSE,
                                      remove_words = "",
                                      second_table = NULL,
                                      max_char_width = 60,
                                      bold_covariates = T)
    }

    return(CCAS_Object)
}
