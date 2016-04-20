#' @title Generate Interpretable Output from a CCAS object.
#' @description Generates a large number of plots of various posterior quantites
#' of interest, including trace plots, parameter estimates, and top words. Note
#' that if the user wishes to generate network plots with nodes colored by some
#' categorical variable, they will need to call the
#' plot_interaction_pattern_network() function separately.
#'
#' @param CCAS_Object The object returned by the ccas() main estimation
#' function.
#' @param output_directory The directory where the user would like to save plots
#' if an output_name_stem is given. Defaults to NULL, in which case the current
#' working directory will be used.
#' @param output_name_stem Defaults to NULL. If not null, then output will be
#' saved to disk using the output_name_stem to differentiate it from output from
#' other model runs.
#' @param generate_plots Logical indicating whether plots should be generated,
#' defaults to TRUE, but may be set to FALSE if the user only wishes to access
#' aggregate data.
#' @return A series of plots and a CCAS_Object containing useful output.
#' @examples
#' \dontrun{
#' # load in saved model output from ccas() function.
#' data(Model_Output)
#' summary_output <- generate_output(Model_Output,
#'                                   output_directory = "~/Desktop",
#'                                   output_name_stem = "Testing")
#' }
#' @export
generate_output <- function(CCAS_Object,
                            output_directory = NULL,
                            output_name_stem = NULL,
                            generate_plots = TRUE) {

    # set the output directory and save the current one so we can reset it
    # afterwards
    if (!is.null(output_directory)) {
        cur_directory <- getwd()
        setwd(output_directory)
    }

    # only generate plots if user specifies that they shoudl be generated
    if (generate_plots) {
        if (!is.null(output_name_stem)) {
            # generate unnormalized topic model log likelihood trace plot and print out
            # Geweke statistic
            pdf(file = paste(output_name_stem,
            "topic_model_log_likelihood.pdf", sep = "_"),
            height = 6,
            width = 8)
            plot_topic_model_log_likelihood(CCAS_Object)
            dev.off()

            # generate LSM trace plots and print out Geweke diagnostics
            pdf(file = paste(output_name_stem,
            "interaction_pattern_log_likelihood.pdf", sep = "_"),
            height = 4,
            width = 8)
            plot_interaction_pattern_log_likelihood(CCAS_Object)
            dev.off()
        } else {
            # generate unnormalized topic model log likelihood trace plot and print out
            # Geweke statistic
            plot_topic_model_log_likelihood(CCAS_Object)
            Sys.sleep(2)

            # generate LSM trace plots and print out Geweke diagnostics
            plot_interaction_pattern_log_likelihood(CCAS_Object)
            Sys.sleep(2)
        }
    }


    # generate interaction pattern specific sub network sociomatrices
    ret_list <- generate_interaction_pattern_subnetworks(
        CCAS_Object)

    # assign output into a list
    CCAS_Object@model_output <- list(
        interaction_pattern_subnetworks = ret_list$interaction_pattern_networks,
        emails_represented = ret_list$emails_represented)

    # generate subnetwork plots for each interaction pattern
    mean_posterior_positions <- vector(mode = "list",
        length = CCAS_Object@interaction_patterns)
    for (i in 1:CCAS_Object@interaction_patterns) {
        if (!is.null(output_name_stem)) {
            pdf(file = paste(output_name_stem,"_interaction_pattern_network_",
                             i,".pdf", sep = ""),
                height = 10,
                width = 10)
            res <- plot_interaction_pattern_network(
                CCAS_Object,
                i,
                plot_color_category = NULL,
                generate_plot = generate_plots)
            dev.off()
        } else {
            res <- plot_interaction_pattern_network(
                CCAS_Object,
                i,
                plot_color_category = NULL,
                generate_plot = generate_plots)
            Sys.sleep(2)
        }
        mean_posterior_positions[[i]] <- res
    }
    # make sure it has interpretable names
    names(mean_posterior_positions) <- paste("interaction_pattern_",
        1:CCAS_Object@interaction_patterns, sep = "")

    # assign the resulting posterior mean positions to out output list object
    CCAS_Object@model_output$mean_posterior_positions <-
        mean_posterior_positions

    # determine a pretty looking number of rows for plot
    ips <- CCAS_Object@interaction_patterns
    if (ips < 7) {
        plots_per_row <- CCAS_Object@interaction_patterns
    } else if (ips == 7 | ips == 8 | ips == 11 | ips == 12) {
        plots_per_row <- 4
    } else {
        plots_per_row <- 5
    }

    # generate plots and get parameter estimates
    pdf(file = paste(output_name_stem,
                     "parameter_estimates.pdf", sep = "_"),
        height = 6,
        width = 12)
    parameter_estimates <- plot_parameter_estimates(CCAS_Object,
                                         plots_per_row = plots_per_row,
                                         generate_plots = generate_plots)
    dev.off()

    CCAS_Object@parameter_estimates <- parameter_estimates

    # get topic and interaction pattern top words
    CCAS_Object <- top_words(CCAS_Object)

    # save the CCAS_Object
    if (!is.null(output_name_stem)) {
        save(CCAS_Object, file = paste(output_name_stem,"model_output.RData",
                                       sep = "_"))
    }

    # reset the working directory
    if (!is.null(output_directory)) {
        setwd(cur_directory)
    }

    return(CCAS_Object)
}
