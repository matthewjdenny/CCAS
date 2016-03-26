generate_output <- function(CCAS_Object,
                            output_directory,
                            output_name_stem) {

    # generate unnormalized topic model log likelihood trace plot and print out
    # Geweke statistic
    plot_topic_model_log_likelihood(CCAS_Object)
    Sys.sleep(2)

    # generate LSM trace plots and print out Geweke diagnostics
    plot_interaction_pattern_log_likelihood(CCAS_Object)
    Sys.sleep(2)

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
        res <- plot_interaction_pattern_network(
            CCAS_Object,
            i,
            plot_color_category = NULL)
        Sys.sleep(2)
        mean_posterior_positions[[i]] <- res
    }
    # make sure it has interpretable names
    names(mean_posterior_positions) <- paste("interaction_pattern_",
        1:CCAS_Object@interaction_patterns, sep = "")

    # assign the resulting posterior mean positions to out output list object
    CCAS_Object@model_output$mean_posterior_positions <-
        mean_posterior_positions

    return(CCAS_Object)
}
