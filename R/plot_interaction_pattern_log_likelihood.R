#' @title Trace plot of interaction pattern model log likelihood
#' @description Generates a trace plot of the interaction pattern model log
#' likelihood. This is a potential way to diagnose the convergence of the model.
#'
#' @param CCAS_Object devtoolAn object of class CCAS containing estimation results.
#' @return A plot.
#' @export
plot_interaction_pattern_log_likelihood <- function(CCAS_Object) {

    # define color
    UMASS_BLUE <- rgb(51,51,153,255,maxColorValue = 255)

    start <- CCAS_Object@final_metropolis_hastings_burnin
    end <- start + CCAS_Object@final_metropolis_hastings_iterations
    prin_seq <- seq(start, end, by = 1/CCAS_Object@thin)[-1]
    # extract the log likelihoods
    log_likelihoods  <- as.numeric(CCAS_Object@MCMC_output$LSM_log_likelihood)

    # set margins
    par(mfrow = c(1,1), mar = c(5,5,4,1))

    # make plot
    plot(y = log_likelihoods,
         x = prin_seq,
         pch = 19,
         col = UMASS_BLUE,
         main = paste("Interaction Pattern Log Likelihood \n",
                      " Geweke Statistic for Last",
                      CCAS_Object@final_metropolis_hastings_iterations,
                      "Iterations:",
                      round(coda::geweke.diag(log_likelihoods)$z,2)),
         xlab = "Iteration",
         ylab = "Log Likelihood",
         cex.lab = 2,
         cex.axis = 1.4,
         cex.main = 1.4)
}
