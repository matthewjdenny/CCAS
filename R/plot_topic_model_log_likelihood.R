#' @title Trace plot of topic model log likelihood
#' @description Generates a trace plot of the unnormalized topic model log
#' likelihood. This is a potential way to diagnose the convergence of the model.
#'
#' @param An object of class CCAS containing estimation results.
#' @param topic_model_burnin The number of iterations of Gibbs sampling to
#' discard before generating trace plot and calculating Geweke statistic.
#' @return A plot.
#' @export
plot_topic_model_log_likelihood <- function(CCAS_Object,
                                            topic_model_burnin = NULL) {

    # define color
    UMASS_BLUE <- rgb(51,51,153,255,maxColorValue = 255)

    # if no burnin is provided, then set it to hald the iterations
    if (is.null(topic_model_burnin)) {
        # if no burnin is provided, remove the first half of iterations
        topic_model_burnin <- floor(CCAS_Object@iterations/2)
    }

    # extract the log likelihoods
    log_likelihoods  <- as.numeric(
        CCAS_Object@topic_model_results$unnoramlized_LDA_log_likelihood)
    len <- length(log_likelihoods)
    geweke_log_likelihoods <- log_likelihoods[topic_model_burnin:len]

    # set margins
    par(mfrow = c(1,1), mar = c(5,5,4,1))

    # make plot
    plot(y = geweke_log_likelihoods,
         x = topic_model_burnin:len,
         pch = 19,
         col = UMASS_BLUE,
         main = paste("Un-Normalized Topic Model Log Likelihood \n",
                      " Geweke Statistic for Last",
                      len - topic_model_burnin,
                      "Iterations:",
                      round(coda::geweke.diag(geweke_log_likelihoods)$z,2)),
         xlab = "Iteration",
         ylab = "Un-Normalized Log Likelihood",
         cex.lab = 2,
         cex.axis = 1.4,
         cex.main = 1.4)
}
