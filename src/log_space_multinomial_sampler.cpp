// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(BH)]]

#include <RcppArmadillo.h>
#include <boost/random.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/detail/disable_warnings.hpp>

using std::log;
using std::exp;
using std::max;

// use the mjd namespace so we can call
namespace mjd {
    int log_space_multinomial_sampler (
            arma::vec unnormalized_discrete_distribution,
            double uniform_draw) {

        // we need to pass in the uniform draw since we only want one RNG in the
        // program and to have it run at the highest level so that we are not
        // constantly re-instantiating the distribution
        // take the log of the uniform draw
        double lud = log(uniform_draw);

        //get the max value of the vector
        double max_val = arma::max(unnormalized_discrete_distribution);

        //get the length of the vector
        int length = unnormalized_discrete_distribution.n_elem;

        //take log sum exp of the entire distribution
        arma::vec exped = exp(unnormalized_discrete_distribution - max_val);
        double lse_dist = max_val + log(arma::sum(exped));

        // now "multiply" (adding in log space) this term by the log uniform
        // draw to get our target value.
        double target = lse_dist + lud;

        // instead of setting to negative infinity, we need to set to zero and
        // then set to the correct value inside the conditional in the for loop
        // so that we do not run into issues with different definitions of the
        // INFINITY constant in different versions of C++.
        double total = 0;

        // int to store the value of what we are going to sample
        int sampled_value = 0;

        //now loop over entries in the distribution
        for (int i = 0; i < length; ++i) {
            // get the value of the current entry in the distribution
            double current = unnormalized_discrete_distribution[i];

            if( i == 0){
                // if we are at the first entry, just set total = current. This
                // works because if total were -Inf, then exp(-Inf - anything) = 0
                // so we would be adding 0 to exp(current - current) = 1, then
                // logging it (= 0), then adding current (= current) so this way
                // we can just skip the step and avoid having to use the negative
                // infinity constant.
                total = current;
            }else{
                // otherwise, we find the current maximum value
                double cur_max = max(current , total);
                // now total becomes log sum exp of those two values
                total = log(exp((current - cur_max)) + exp((total - cur_max))) + cur_max;
            }

            //if our total is now greater than the target, then we have found the
            //right index so save it, break the loop, and return the value
            if (total > target){
                sampled_value = i;
                break;
            }
        }
        return(sampled_value);
    }
}

// for testing we will wrap and export this function so it is available in R.
// [[Rcpp::export]]
int lsms(arma::vec unnormalized_discrete_distribution,
         int seed){

    // Set RNG and define uniform distribution
    boost::mt19937 generator(seed);
    boost::uniform_01<double> uniform_distribution;

    // get the random uniform draw and log it
    double rand_num = uniform_distribution(generator);

    // take a draw from the unnormalized log distribution
    int temp = mjd::log_space_multinomial_sampler (
            unnormalized_discrete_distribution,
            rand_num);
    //increment by one to go back to 1 indexing in R.
    return(temp + 1);
}
