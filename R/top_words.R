#' @title Generate Topic top words
#' @description Generates data.frames containing the top words associated with
#' each topic and with each interaction pattern
#'
#' @param CCAS_Object The object returned by the ccas() main estimation
#' function.
#' @param output_LaTeX_code Logical indicating whether LaTeX code should be
#' output by the functon for incorporating in tables.
#' @return A list of data.frames with topic and interaction pattern top words.
#' @export
top_words <- function(CCAS_Object,
                      output_LaTeX_code) {

    # allocate some variables we will use below:
    interaction_patterns <- CCAS_Object@interaction_patterns
    topics <- CCAS_Object@number_of_topics
    vocab_size <- CCAS_Object@ComNet_Object@vocabulary_size
    vocab <- CCAS_Object@ComNet_Object@vocabulary
    topic_ips <- CCAS_Object@topic_model_results$topic_interaction_patterns
    topic_ips <- topic_ips[nrow(topic_ips),]
    topic_word_counts <- CCAS_Object@topic_model_results$word_type_topic_counts

    # allocate a dataframe for top words and another for ordered counts
    topic_top_words <- matrix(NA, nrow = topics, ncol = vocab_size)
    topic_top_word_counts <- matrix(0, nrow = topics, ncol = vocab_size)

    # loop through topics and generate ordered counts
    for (i in 1:topics) {

    }




}
