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
    Rcpp::List sample_new_interaction_pattern_parameters(
        arma::vec intercepts,
        arma::mat coefficients,
        arma::cube latent_positions,
        double intercept_prior_variance,
        double coefficient_prior_variance,
        double latent_position_prior_variance,
        bool using_coefficients) {

        // get number of interaction patterns
        int num_interaction_patterns = intercepts.n_elem;

        //get the number of if we are using them, set this eqaul to two since we
        //need to allocate the matrix even if we are not using it.
        int num_coefficients = 2;
        if(using_coefficients){
            num_coefficients = coefficients.n_cols;
        }

        //get the number of actors
        int num_actors = latent_positions.n_cols;

        //get the number of latent dimensions
        int num_latent_dimensions = latent_positions.n_slices;

        // initialize the proposed values
        arma::vec proposed_intercepts = arma::zeros(num_interaction_patterns);
        arma::mat proposed_coefficients = arma::zeros(num_interaction_patterns,
                                                      num_coefficients);
        arma::cube proposed_latent_positions = arma::zeros(
            num_interaction_patterns,
            num_actors,
            num_latent_dimensions);

        //sample new intercepts centered at current ones.
        for (int i = 0; i < num_interaction_patterns; ++i) {
            // see https://svn.r-project.org/R/trunk/src/nmath/rnorm.c for source
            // code. First argument is the mean of the normal distribution we
            // are sampling from and second is its variance.
            proposed_intercepts[i]= R::rnorm(intercepts[i],
                                              intercept_prior_variance);
        }

        // if we are using coefficients, sample new coefficients centered at
        // current ones.
        if (using_coefficients) {
            for (int i = 0; i < num_interaction_patterns; ++i) {
                for (int j = 0; j < num_coefficients; ++j) {
                    // see https://svn.r-project.org/R/trunk/src/nmath/rnorm.c for source
                    // code. First argument is the mean of the normal distribution we
                    // are sampling from and second is its variance.
                    proposed_coefficients(i,j) += R::rnorm(coefficients(i,j),
                                                    coefficient_prior_variance);
                }
            }
        }

        //sample new latent positions centered at current ones.
        for (int i = 0; i < num_interaction_patterns; ++i) {
            for (int j = 0; j < num_actors; ++j) {
                for (int k = 0; k < num_latent_dimensions; ++k) {
                    // see https://svn.r-project.org/R/trunk/src/nmath/rnorm.c for source
                    // code. First argument is the mean of the normal distribution we
                    // are sampling from and second is its variance.
                    proposed_latent_positions(i,j,k) = R::rnorm(
                        latent_positions(i,j,k),
                        latent_position_prior_variance);
                }
            }
        }

        Rcpp::List to_return(3);
        to_return[0] = proposed_intercepts;
        to_return[1] = proposed_coefficients;
        to_return[2] = proposed_latent_positions;
        return to_return;
    }
}

using namespace Rcpp;
// for testing we will wrap and export this function so it is available in R.
// [[Rcpp::export]]
List snipp(arma::vec intercepts,
          arma::mat coefficients,
          NumericVector latent_pos,
          double intercept_prior_variance,
          double coefficient_prior_variance,
          double latent_position_prior_variance,
          bool using_coefficients){

    // we have to do this stupid trick to pass in 3d arrays from R. We pass in as
    // a vector, then instatiate a cube object from there.
    IntegerVector arrayDims = latent_pos.attr("dim");
    arma::cube latent_positions(latent_pos.begin(), arrayDims[0], arrayDims[1],
                                arrayDims[2], false);

    List new_params = mjd::sample_new_interaction_pattern_parameters (
        intercepts,
        coefficients,
        latent_positions,
        intercept_prior_variance,
        coefficient_prior_variance,
        latent_position_prior_variance,
        using_coefficients);

    return new_params;
}

