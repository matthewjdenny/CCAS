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




}
