// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(BH)]]

#include <RcppArmadillo.h>
#include <boost/random.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/detail/disable_warnings.hpp>

using std::log;
using std::exp;
using std::max;
using std::abs;
// use the mjd namespace so we can call in other functions
namespace mjd {
    double prior_probability_interaction_pattern_parameters(
            arma::vec intercepts,
            arma::mat coefficients,
            arma::cube latent_positions,
            double intercept_prior_mean,
            double intercept_prior_variance,
            double coefficient_prior_mean,
            double coefficient_prior_variance,
            double latent_position_prior_mean,
            double latent_position_prior_variance,
            bool using_coefficients) {

        // get number of interaction patterns
        int num_interaction_patterns = intercepts.n_elem;

        //get the number of if we are using them
        int num_coefficients = 0;
        if(using_coefficients){
            num_coefficients = coefficients.n_cols;
        }

        //get the number of actors
        int num_actors = latent_positions.n_cols;

        //get the number of latent dimensions
        int num_latent_dimensions = latent_positions.n_slices;

        // initialize the prior probability to zero
        double log_prior_probability = 0;


        //first get all of the intercept prior probabilities
        for (int i = 0; i < num_interaction_patterns; ++i) {
            // see https://svn.r-project.org/R/trunk/src/nmath/pnorm.c for source
            // code. First argument is the value we want to evaluate, second is prior
            // mean, third is prior variance, fourth is whether function should
            // return the log probability (1), which is what we want.
            log_prior_probability += R::dnorm(intercepts[i],
                                              intercept_prior_mean,
                                              intercept_prior_variance,
                                              1);
        }

        // if we are using coefficients, then get their prior probabilities.
        if (using_coefficients) {
            for (int i = 0; i < num_interaction_patterns; ++i) {
                for (int j = 0; j < num_coefficients; ++j) {
                    // see https://svn.r-project.org/R/trunk/src/nmath/pnorm.c for
                    // source code. First argument is the value we want to evaluate,
                    // second is prior mean, third is prior variance, fourth is
                    // whether function should return the log probability (1), which
                    // is what we want.
                    log_prior_probability += R::dnorm(coefficients(i,j),
                                                      coefficient_prior_mean,
                                                      coefficient_prior_variance,
                                                      1);
                }
            }
        }

        // now for the latent positions
        for (int i = 0; i < num_interaction_patterns; ++i) {
            for (int j = 0; j < num_actors; ++j) {
                for (int k = 0; k < num_latent_dimensions; ++k) {
                    // see https://svn.r-project.org/R/trunk/src/nmath/pnorm.c for
                    // source code. First argument is the value we want to evaluate,
                    // second is prior mean, third is prior variance, fourth is
                    // whether function should return the log probability (1), which
                    // is what we want.
                    log_prior_probability += R::dnorm(latent_positions(i,j,k),
                                                      latent_position_prior_mean,
                                                      latent_position_prior_variance,
                                                      1);
                }
            }
        }

        return log_prior_probability;
    }
}

using namespace Rcpp;
// for testing we will wrap and export this function so it is available in R.
// [[Rcpp::export]]
int ppipp(arma::vec intercepts,
          arma::mat coefficients,
          NumericVector latent_pos,
          double intercept_prior_mean,
          double intercept_prior_variance,
          double coefficient_prior_mean,
          double coefficient_prior_variance,
          double latent_position_prior_mean,
          double latent_position_prior_variance,
          bool using_coefficients){

    // we have to do this stupid trick to pass in 3d arrays from R. We pass in as
    // a vector, then instatiate a cube object from there.
    IntegerVector arrayDims = latent_pos.attr("dim");
    arma::cube latent_positions(latent_pos.begin(), arrayDims[0], arrayDims[1],
                                  arrayDims[2], false);

    double log_prob = mjd::prior_probability_interaction_pattern_parameters (
        intercepts,
        coefficients,
        latent_positions,
        intercept_prior_mean,
        intercept_prior_variance,
        coefficient_prior_mean,
        coefficient_prior_variance,
        latent_position_prior_mean,
        latent_position_prior_variance,
        using_coefficients);

    return log_prob;
}
