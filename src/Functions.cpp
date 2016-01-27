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
        return sampled_value;
    }
	
	
    double sum_over_t_edge_probability (
            arma::cube edge_probabilities,
            int tokens_in_document,
            int current_token_topic_assignment,
            arma::vec current_document_topic_counts,
            bool leave_out_current_token,
            arma::vec topic_interaction_patterns,
            int document_sender,
            int document_recipient,
            int leave_out_topic) {

        // get number of topics
        int num_topics = current_document_topic_counts.n_elem;

        // initialize sum to zero
        double sum = 0;

        // if we are leaving out the current token then decrement its count.
        if (leave_out_current_token) {
            current_document_topic_counts[current_token_topic_assignment] -= 1;
        }

        //  leave out topic can default to -1 in which case we do not leave out
        //  any topics. This outer conditional makes things run faster for most
        //  cases.
        if (leave_out_topic > -1) {
            for (int i = 0; i < num_topics; ++i) {
                if (leave_out_topic != i) {
                    sum += (double(current_document_topic_counts[i])/
                                double(tokens_in_document)) *
                                    edge_probabilities(document_sender,
                                                       document_recipient,
                                                       i);
                }
            }
        } else {
            for (int i = 0; i < num_topics; ++i) {
                sum += (double(current_document_topic_counts[i])/
                            double(tokens_in_document)) *
                                edge_probabilities(document_sender,
                                                   document_recipient,
                                                   i);
            }
        }

        return sum;
    }

    double lsm_contribution(
            arma::cube edge_probabilities,
            int tokens_in_document,
            int topic,
            int current_token_topic_assignment,
            arma::vec current_document_topic_counts,
            arma::vec document_edge_values,
            arma::vec topic_interaction_patterns,
            int document_sender,
            int current_topic) {

        // get number of interaction patterns
        int num_actors = document_edge_values.n_elem;

        // initialize the lsm contribution
        double contribution = 0;

        //first get all of the intercept prior probabilities
        for (int i = 0; i < num_actors; ++i) {
            if (i != document_sender) {
                double sum_term = sum_over_t_edge_probability (
                    edge_probabilities,
                    tokens_in_document,
                    current_token_topic_assignment,
                    current_document_topic_counts,
                    -1,
                    topic_interaction_patterns,
                    document_sender,
                    i,
                    false);
                //make sure we do not do integer division so we cast as a double
                sum_term += 1/(double(tokens_in_document)) *
                    edge_probabilities(document_sender,
                                       i,
                                       topic_interaction_patterns[topic]);

                //now we take the log of the sum term or 1- sum term depending
                //on whether the edge value is equal to 1 (present), or 0
                //(absent).
                if (document_edge_values[i] == 1) {
                    contribution += log(sum_term);
                } else {
                    // we can use this equality to simplify things:
                    // (1-a)+(1-b) = 2-a-b = 2-(a+b)
                    contribution += log(2 - sum_term);
                }
            }
        }
        return contribution;
    }
	
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
double ep(arma::vec intercepts,
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
    return (temp + 1);
}

// [[Rcpp::export]]
int sotep(NumericVector edge_probs,
          int tokens_in_document,
          int current_token_topic_assignment,
          arma::vec current_document_topic_counts,
          bool leave_out_current_token,
          arma::vec topic_interaction_patterns,
          int document_sender,
          int document_recipient,
          int leave_out_topic){

    // we have to do this stupid trick to pass in 3d arrays from R. We pass in as
    // a vector, then instatiate a cube object from there.
    IntegerVector arrayDims = edge_probs.attr("dim");
    arma::cube edge_probabilities(edge_probs.begin(), arrayDims[0], arrayDims[1],
                                arrayDims[2], false);

    double sum = mjd::sum_over_t_edge_probability (
        edge_probabilities,
        tokens_in_document,
        current_token_topic_assignment,
        current_document_topic_counts,
        leave_out_current_token,
        topic_interaction_patterns,
        document_sender,
        document_recipient,
        leave_out_topic);

    return sum;
}

// [[Rcpp::export]]
double lsmc(NumericVector edge_probs,
             int tokens_in_document,
             int topic,
             int current_token_topic_assignment,
             arma::vec current_document_topic_counts,
             arma::vec document_edge_values,
             arma::vec topic_interaction_patterns,
             int document_sender,
             int current_topic){

    // we have to do this stupid trick to pass in 3d arrays from R. We pass in as
    // a vector, then instatiate a cube object from there.
    IntegerVector arrayDims = edge_probs.attr("dim");
    arma::cube edge_probabilities(edge_probs.begin(), arrayDims[0], arrayDims[1],
                                arrayDims[2], false);

    double contrib = mjd::lsm_contribution (
        edge_probabilities,
        tokens_in_document,
        topic,
        current_token_topic_assignment,
        current_document_topic_counts,
        document_edge_values,
        topic_interaction_patterns,
        document_sender,
        current_topic);

    return contrib;
}

// [[Rcpp::export]]
double ppipp(arma::vec intercepts,
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
