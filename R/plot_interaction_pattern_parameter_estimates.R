#' @title Plot interaction-pattern sepcific coefficients.
#' @description Plots parameter estimates and intercept term associated with a
#' particular interaction pattern.
#'
#' @param CCAS_Object The object returned by the ccas() main estimation
#' function.
#' @param interaction_pattern_index The index of the cluster specific
#' communication network we wish to plot.
#' @param coefficient_names Defualts to NULL. Can be a string vector of names
#' for coefficients to be used in making publication quality plots.
#' @param leave_out_coefficients Defaults to NULL. Can be a string vector of
#' coefficient names as they appear in the plot. These coefficients will be
#' removed from the final plot. Useful if the intercept term is much larger in
#' magnitude than other estimates, and the user wishes to clarify the other
#' parameter estimates without normalizing.
#' @param normalize_coefficients Defaults to FALSE, if TRUE then parameter
#' estimates will be converted be deivided by their standard deviations with
#' and displayed with 95 percent confidence intervals. These coefficients will
#' no longer be comparable, but make graphical interpretation of significance
#' and sign easier.
#' @param generate_plot Logical indicating whether a plot should be generated,
#' defaults to TRUE, but may be set to FALSE if the user only wishes to access
#' parameter estimates.
#' @return A plot.
#' @export
plot_interaction_pattern_parameter_estimates <- function(
    CCAS_Object,
    interaction_pattern_index,
    coefficient_names = NULL,
    leave_out_coefficients = NULL,
    normalize_coefficients = FALSE,
    generate_plot = TRUE) {

    UMASS_BLUE <- rgb(51,51,153,255,maxColorValue = 255)

    # get the number of samples after burnin and thinning
    intercept_mean <- mean(CCAS_Object@MCMC_output$intercepts[
        ,interaction_pattern_index])
    intercept_sd <- sd(CCAS_Object@MCMC_output$intercepts[
        ,interaction_pattern_index])

    # determine whether we are using covariates
    using_covariates <- CCAS_Object@number_of_covariates > 0

    # if we are using covariates, get the mean and sd of coefficient
    # parameter estimates and add these to the model data frame
    if (using_covariates) {
        coefficients <- CCAS_Object@MCMC_output$coefficients[
            interaction_pattern_index,,]
        coefficient_means <- apply(coefficients,1,mean)
        coefficient_sds <- apply(coefficients,1,sd)
        modelFrame <- data.frame(
            Variable = c("intercept",
                         dimnames(CCAS_Object@covariate_array)[[3]]) ,
            Coefficient = c(intercept_mean,
                            coefficient_means),
            SE = c(intercept_sd,
                   coefficient_sds),
            Interaction_Pattern = rep(
                paste("Interaction Pattern",interaction_pattern_index),
                length(c(intercept_sd,coefficient_sds))),
            Color = rep("Blue",
            length(c(intercept_sd,coefficient_sds))),
            stringsAsFactors = FALSE)
        data <- data.frame(modelFrame)
    } else {
        # only plot the intercept
        modelFrame <- data.frame(
            Variable = "intercept",
            Coefficient = intercept_mean,
            SE = intercept_sd,
            Interaction_Pattern = paste("Interaction Pattern",
                                        interaction_pattern_index),
            Color = "Blue",
            stringsAsFactors = FALSE)
        data <- data.frame(modelFrame)
    }


    # rename coefficients if necessary
    if (!is.null(coefficient_names)) {
        if (length(coefficient_names) != nrow(data)) {
            stop("coefficient_names must be the same length as the number of covariates in the plot. Try setting coefficient_names = NULL and counting the number of coefficients to check that you have provided the right number.")
        }
        cat("Replacing:\n")
        print(data$Variable)
        cat("With..\n")
        print(coefficient_names)
        data$Variable <- coefficient_names
    }

    #remove variables
    if (!is.null(leave_out_coefficients)) {
        for (j in 1:length(leave_out_coefficients)) {
            remove <- which(data$Variable == leave_out_coefficients[j])
            if (length(remove) == 1) {
                data <- data[-remove,]
                cat("Removing variable",leave_out_coefficients[j],"\n")
            } else if (length(remove) > 1) {
                cat("Your argument matches more than one variable. Please respecify.\n")
            } else {
                cat ("There is no variable",leave_out_coefficients[j],
                     ".Please respecify.\n")
            }
        }
    }

    # standardize coefficients
    if (normalize_coefficients) {
        for (i in 1:nrow(data)) {
            data$Coefficient[i] <- data$Coefficient[i]/data$SE[i]
            data$SE[i] <- 1
        }
    }


    # generate plot
    zp1 <- ggplot2::ggplot(data, ggplot2::aes(colour = Color)) +
           ggplot2::scale_color_manual(values = UMASS_BLUE)

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

    if (generate_plot) {
        print(zp1)
    }

    # return everything
    return(data)
}
