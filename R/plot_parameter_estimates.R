#' @title Plot interaction-pattern sepcific sub-networks
#' @description Plots coefficient ad intercept estimates for all interaction
#' patterns.
#'
#' @param CCAS_Object The object returned by the ccas() main estimation
#' function.
#' @param plots_per_row The number of interaction pattern coefficient plots to
#' include in each row of the output figure.
#' @param generate_plots Logical indicating whether plots should be generated,
#' defaults to TRUE, but may be set to FALSE if the user only wishes to access
#' parameter estimates.
#' @return A plot.
#' @examples
#' \dontrun{
#' # load in saved model output from ccas() function.
#' data(Model_Output)
#' # plots_per_row should be equal to the number of interactions patterns.
#' summary_data <- plot_parameter_estimates(Model_Output, plots_per_row = 4)
#' }
#' @export
plot_parameter_estimates <- function(CCAS_Object,
                                     plots_per_row,
                                     generate_plots = TRUE){

    Interaction_Pattern <- Color <- Variable <- Coefficient <- SE <- NULL
    # loop over interaction patterns
    for (i in  1:CCAS_Object@interaction_patterns) {
        temp <- plot_interaction_pattern_parameter_estimates(
            CCAS_Object,
            interaction_pattern_index = i,
            coefficient_names = NULL,
            leave_out_coefficients = NULL,
            normalize_coefficients = FALSE,
            generate_plot = FALSE)
        if (i == 1) {
            data <- temp
        } else {
            data <- rbind(data,temp)
        }

    }

    # generate plot
    UMASS_BLUE <- rgb(51,51,153,255,maxColorValue = 255)
    zp1 <- ggplot2::ggplot(data, ggplot2::aes(colour = Color)) +
        ggplot2::scale_color_manual(values = UMASS_BLUE) +
        ggplot2::facet_wrap(~Interaction_Pattern, ncol = plots_per_row)
    zp1 <- zp1 + ggplot2::geom_hline(
        yintercept = 0,
        colour = gray(1/2),
        lty = 2)
    zp1 <- zp1 + ggplot2::geom_linerange( ggplot2::aes(
        x = Variable,
        ymin = Coefficient - SE*(-qnorm((1 - 0.9)/2)),
        ymax = Coefficient + SE*(-qnorm((1 - 0.9)/2))),
        lwd = 1,
        position = ggplot2::position_dodge(width = 1/2))
    zp1 <- zp1 + ggplot2::geom_pointrange(ggplot2::aes(
        x = Variable,
        y = Coefficient,
        ymin = Coefficient - SE*(-qnorm((1 - 0.95)/2)),
        ymax = Coefficient + SE*(-qnorm((1 - 0.95)/2))),
        lwd = 1/2,
        position = ggplot2::position_dodge(width = 1/2),
        shape = 21, fill = "WHITE")

    zp1 <- zp1  + ggplot2::theme_bw() +
        ggplot2::coord_flip() +
        ggplot2::theme(legend.position = "none")

    if (generate_plots) {
        print(zp1)
    }

    # remove the coloring column
    data <- data[,-5]

    return(data)
}
