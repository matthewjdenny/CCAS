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
    double edge_probability (
            arma::vec intercepts,
            arma::mat coefficients,
            arma::cube latent_positions,
            int sender,
            int recipient,
            arma::vec current_covariates,
            int interaction_pattern,
            bool using_coefficients) {

        // initialize linear predictor using intercept
        double linear_predictor = intercepts[interaction_pattern];

        // if we are using coefficients, then go and add their effects
        if (using_coefficients) {
            //now loop over coefficients for current interaction pattern
            int num_coefs = coefficients.n_cols;
            for (int i = 0; i < num_coefs; ++i) {
                linear_predictor += coefficients(interaction_pattern,i) *
                    current_covariates[i];
            }
        }

        // subtract out distance between two actors in latent space
        int num_latent_spaces = latent_positions.n_slices;
        for (int i = 0; i < num_latent_spaces; ++i) {
            linear_predictor += abs((latent_positions(interaction_pattern, sender, i) -
                latent_positions(interaction_pattern, recipient, i)));
        }

        // now run through logit
        double edge_prob = (1/double(1+ exp(-linear_predictor)));
        return edge_prob;
    }
}

using namespace Rcpp;
// for testing we will wrap and export this function so it is available in R.
// [[Rcpp::export]]
int ep(arma::vec intercepts,
       arma::mat coefficients,
       NumericVector latent_pos,
       int sender,
       int recipient,
       arma::vec current_covariates,
       int interaction_pattern,
       bool using_coefficients){

    // we have to do this stupid trick to pass in 3d arrays from R. We pass in as
    // a vector, then instatiate a cube object from there.
    IntegerVector arrayDims = latent_pos.attr("dim");
    arma::cube latent_positions(latent_pos.begin(), arrayDims[0], arrayDims[1],
                                arrayDims[2], false);

    double prob = mjd::edge_probability (intercepts,
            coefficients,
            latent_positions,
            sender,
            recipient,
            current_covariates,
            interaction_pattern,
            using_coefficients);

    return prob;
}
