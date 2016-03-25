generate_output <- function(CCAS_Object,
                            output_directory,
                            output_name_stem) {

    # generate unnormalized topic model log likelihood trace plot and print out
    # Geweke statistic
    plot_topic_model_log_likelihood(CCAS_Object)

    # generate LSM trace plots and print out Geweke diagnostics
    plot_interaction_pattern_log_likelihood(CCAS_Object)

    # generate interaction pattern specific sub network sociomatrices
    interaction_pattern_subnetworks <- generate_interaction_pattern_subnetworks(
        CCAS_Object)

    # assign output into a list
    CCAS_Object@model_output <- list(
        interaction_pattern_subnetworks = interaction_pattern_subnetworks)
    return(CCAS_Object)
}
