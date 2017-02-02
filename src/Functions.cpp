#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(BH)]]
// [[Rcpp::depends(RcppParallel)]]

#include <RcppParallel.h>
#include <RcppArmadillo.h>
#include <boost/random.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/detail/disable_warnings.hpp>

using std::log;
using std::exp;
using std::max;
using std::abs;
using std::sqrt;
using std::pow;
// use the mjd namespace so we can call in other functions
namespace mjd {
    using namespace Rcpp ;
    using std::log;
    using std::exp;
    using std::max;
    using std::abs;
    using std::sqrt;
    using std::pow;

    // ***********************************************************************//
    //                 Log Space Multinomial Sampler                          //
    // ***********************************************************************//

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

    // ***********************************************************************//
    //                       Edge Probability                                 //
    // ***********************************************************************//

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
        double squared_dist = 0;
        //generate the sum of squares
        for (int i = 0; i < num_latent_spaces; ++i) {
            squared_dist += pow(
                (latent_positions(interaction_pattern, sender, i) -
                latent_positions(interaction_pattern, recipient, i)), 2);
        }
        // take the square root and subtract it (the further away nodes are, the
        // less likely the tie.)
        linear_predictor -= sqrt(squared_dist);

        // for testing, print out the linear predictor
        //Rcpp::Rcout << "Linear Predictor: " << linear_predictor << std::endl;

        // now run through logit
        double edge_prob = (1/double(1 + exp(-linear_predictor)));
        return edge_prob;
    }

	// ***********************************************************************//
	//                 Sum Over T Edge Probabilities                          //
	// ***********************************************************************//

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
                                    edge_probabilities(
                                        document_sender,
                                        document_recipient,
                                        topic_interaction_patterns[i]);
                }
            }
        } else {
            for (int i = 0; i < num_topics; ++i) {
                sum += (double(current_document_topic_counts[i])/
                            double(tokens_in_document)) *
                                edge_probabilities(
                                    document_sender,
                                    document_recipient,
                                    topic_interaction_patterns[i]);
                // for testing, print out the sum as we generate it.
                // Rcpp::Rcout << "Sum: " << sum << std::endl;
            }
        }

        return sum;
    }

	// ***********************************************************************//
	//         Get Prior Probability of Interaction Pattern Parameters        //
	// ***********************************************************************//

    double prior_probability_interaction_pattern_parameters(
            arma::vec intercepts,
            arma::mat coefficients,
            arma::cube latent_positions,
            double intercept_prior_mean,
            double intercept_prior_standard_deviation,
            double coefficient_prior_mean,
            double coefficient_prior_standard_deviation,
            double latent_position_prior_mean,
            double latent_position_prior_standard_deviation,
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
            // mean, third is prior standard_deviation, fourth is whether function should
            // return the log probability (1), which is what we want.
            log_prior_probability += R::dnorm(intercepts[i],
                                              intercept_prior_mean,
                                              intercept_prior_standard_deviation,
                                              1);
        }

        // if we are using coefficients, then get their prior probabilities.
        if (using_coefficients) {
            for (int i = 0; i < num_interaction_patterns; ++i) {
                for (int j = 0; j < num_coefficients; ++j) {
                    // see https://svn.r-project.org/R/trunk/src/nmath/pnorm.c for
                    // source code. First argument is the value we want to evaluate,
                    // second is prior mean, third is prior standard_deviation, fourth is
                    // whether function should return the log probability (1), which
                    // is what we want.
                    log_prior_probability += R::dnorm(coefficients(i,j),
                                                      coefficient_prior_mean,
                                                      coefficient_prior_standard_deviation,
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
                    // second is prior mean, third is prior standard_deviation, fourth is
                    // whether function should return the log probability (1), which
                    // is what we want.
                    log_prior_probability += R::dnorm(latent_positions(i,j,k),
                                                      latent_position_prior_mean,
                                                      latent_position_prior_standard_deviation,
                                                      1);
                }
            }
        }

        return log_prior_probability;
    }

    // ***********************************************************************//
    //                Sample New Interaction Pattern Parameters               //
    // ***********************************************************************//

    Rcpp::List sample_new_interaction_pattern_parameters(
        arma::vec intercepts,
        arma::mat coefficients,
        arma::cube latent_positions,
        arma::vec intercept_proposal_standard_deviations,
        arma::vec coefficient_proposal_standard_deviations,
        arma::vec latent_position_proposal_standard_deviations,
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
            // are sampling from and second is its standard_deviation.
            proposed_intercepts[i] = R::rnorm(intercepts[i],
                                              intercept_proposal_standard_deviations[i]);
        }

        // if we are using coefficients, sample new coefficients centered at
        // current ones.
        if (using_coefficients) {
            for (int i = 0; i < num_interaction_patterns; ++i) {
                for (int j = 0; j < num_coefficients; ++j) {
                    // see https://svn.r-project.org/R/trunk/src/nmath/rnorm.c for source
                    // code. First argument is the mean of the normal distribution we
                    // are sampling from and second is its standard_deviation.
                    proposed_coefficients(i,j) = R::rnorm(coefficients(i,j),
                                                    coefficient_proposal_standard_deviations[i]);
                }
            }
        }

        //sample new latent positions centered at current ones.
        for (int i = 0; i < num_interaction_patterns; ++i) {
            for (int j = 0; j < num_actors; ++j) {
                for (int k = 0; k < num_latent_dimensions; ++k) {
                    // see https://svn.r-project.org/R/trunk/src/nmath/rnorm.c for source
                    // code. First argument is the mean of the normal distribution we
                    // are sampling from and second is its standard_deviation.
                    proposed_latent_positions(i,j,k) = R::rnorm(
                        latent_positions(i,j,k),
                        latent_position_proposal_standard_deviations[i]);
                }
            }
        }

        Rcpp::List to_return(3);
        to_return[0] = proposed_intercepts;
        to_return[1] = proposed_coefficients;
        to_return[2] = proposed_latent_positions;
        return to_return;
    }

    // ***********************************************************************//
    //                           LSM Contribution                             //
    // ***********************************************************************//

    double lsm_contribution(
            arma::cube edge_probabilities,
            int tokens_in_document,
            int topic,
            int current_token_topic_assignment,
            arma::vec current_document_topic_counts,
            arma::vec document_edge_values,
            arma::vec topic_interaction_patterns,
            int document_sender) {

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
                    true,
                    topic_interaction_patterns,
                    document_sender,
                    i,
                    -1);
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


    // ***********************************************************************//
    //                           LDA Contribution                             //
    // ***********************************************************************//

    double lda_contribution(
            int tokens_in_document,
            int current_token_topic_assignment,
            arma::vec current_document_topic_counts,
            arma::mat word_type_topic_counts,
            arma::vec topic_token_counts,
            int topic,
            int current_word_type,
            arma::vec alpha_m,
            arma::vec beta_n,
            double beta) {

        int current_dtc = current_document_topic_counts[topic];

        // get the word-type topic count for the current word and topic
        int wttc = word_type_topic_counts(current_word_type,topic);

        // get the total count of tokens currently associated with topic
        int tc = topic_token_counts[topic];

        // the topic is the same as the current assignement, then we decrement
        // the relevant counts since we are doing Gibbs sampling so we need to
        // hold out hte current position
        if (topic == current_token_topic_assignment) {
            current_dtc -= 1;
            wttc -= 1;
            tc -= 1;
        }

        // now we calculate the contribution (which will be in log space)
        double contribution = log(double(current_dtc)  + alpha_m[topic]) +
            log(wttc + beta_n[current_word_type]) - log(tc + beta);

        // return the value
        return contribution;
    }

    // ***********************************************************************//
    //     Constructor for Parallel Token Topic Distribution Generator        //
    // ***********************************************************************//

    // create a RcppParallel::Worker struct that we can use to fill in the
    // entries in our token topic distribution in parallel
    struct Parallel_Token_Topic_Distribution : public RcppParallel::Worker {

        // Instantiate all of our input variables which will then be initialized
        // via the Function initializer below. This is how the struct makes sure
        // that everything has the right type
        arma::cube edge_probabilities;
        int tokens_in_document;
        int current_token_topic_assignment;
        arma::vec current_document_topic_counts;
        arma::mat word_type_topic_counts;
        arma::vec topic_token_counts;
        int current_word_type;
        arma::vec alpha_m;
        arma::vec beta_n;
        arma::vec document_edge_values;
        arma::vec topic_interaction_patterns;
        int document_sender;
        double beta;

        // We need to initialize the output vector to an RcppParallel::RVector
        // vector (which is compatible with an Rcpp::NumericVector but not an
        // arma::vec). This is what will be implicitly returned
        RcppParallel::RVector<double> return_dist;

        // Function initializer
        Parallel_Token_Topic_Distribution(arma::cube edge_probabilities,
                    int tokens_in_document,
                    int current_token_topic_assignment,
                    arma::vec current_document_topic_counts,
                    arma::mat word_type_topic_counts,
                    arma::vec topic_token_counts,
                    int current_word_type,
                    arma::vec alpha_m,
                    arma::vec beta_n,
                    arma::vec document_edge_values,
                    arma::vec topic_interaction_patterns,
                    int document_sender,
                    double beta,
                    Rcpp::NumericVector return_dist)
        : edge_probabilities(edge_probabilities),
        tokens_in_document(tokens_in_document),
        current_token_topic_assignment(current_token_topic_assignment),
        current_document_topic_counts(current_document_topic_counts),
        word_type_topic_counts(word_type_topic_counts),
        topic_token_counts(topic_token_counts),
        current_word_type(current_word_type),
        alpha_m(alpha_m),
        beta_n(beta_n),
        document_edge_values(document_edge_values),
        topic_interaction_patterns(topic_interaction_patterns),
        document_sender(document_sender),
        beta(beta),
        return_dist(return_dist) {}

        // function call operator that works for the specified range (begin/end)
        void operator()(std::size_t begin, std::size_t end) {
            for (std::size_t i = begin; i < end; i++) {
                // this is the exact same code as is in the non-parallel version
                // (I copy-pasted it).
                double lsm_contr = mjd::lsm_contribution (
                    edge_probabilities,
                    tokens_in_document,
                    i,
                    current_token_topic_assignment,
                    current_document_topic_counts,
                    document_edge_values,
                    topic_interaction_patterns,
                    document_sender);
                // for testing.
                // Rcpp::Rcout << "LSM Contribution: " << lsm_contr << std::endl;

                double lda_contr = mjd::lda_contribution(
                    tokens_in_document,
                    current_token_topic_assignment,
                    current_document_topic_counts,
                    word_type_topic_counts,
                    topic_token_counts,
                    i,
                    current_word_type,
                    alpha_m,
                    beta_n,
                    beta);
                // for testing.
                // Rcpp::Rcout << "LDA Contribution: " << lda_contr << std::endl;

                // add the values (since we are working in log space) and put them
                // in the appropriate bin the in the distribution.
                return_dist[i] = lsm_contr + lda_contr;
            }
        }
    };

    // ***********************************************************************//
    //             Calculate token topic probabilities in parallel            //
    // ***********************************************************************//
    arma::vec parallel_token_topic_probabilities(
            arma::cube edge_probabilities,
            int tokens_in_document,
            int current_token_topic_assignment,
            arma::vec current_document_topic_counts,
            arma::mat word_type_topic_counts,
            arma::vec topic_token_counts,
            int current_word_type,
            arma::vec alpha_m,
            arma::vec beta_n,
            arma::vec document_edge_values,
            arma::vec topic_interaction_patterns,
            int document_sender,
            double beta,
            int number_of_topics) {

        // the vector that will be operated on by RcppParallel::parallelFor.
        // This vector must be an Rcpp::NumericVector or the parallelization
        // will not work. THis makes everything slower, but I have not found an
        // alternative
        Rcpp::NumericVector output_vec(number_of_topics);
        // Once we have the token topic distribution back from the parallel for
        // function, we will put it in to an arma::vec of the same length.
        arma::vec return_vec = arma::zeros(number_of_topics);

        // create the RcppParallel::Worker function which will iterate over the
        // distribution
        Parallel_Token_Topic_Distribution Parallel_Token_Topic_Distribution(
                edge_probabilities,
                tokens_in_document,
                current_token_topic_assignment,
                current_document_topic_counts,
                word_type_topic_counts,
                topic_token_counts,
                current_word_type,
                alpha_m,
                beta_n,
                document_edge_values,
                topic_interaction_patterns,
                document_sender,
                beta,
                output_vec);

        // Call our Parallel_Token_Topic_Distribution with parallelFor
        RcppParallel::parallelFor(0,
                                  number_of_topics,
                                  Parallel_Token_Topic_Distribution);

        // Take the NumericVector output and put it in the arma::vec to actually
        // return
        for (int i = 0; i < number_of_topics; i++) {
            // write to output vector
            return_vec[i] = output_vec[i];
        }

        return return_vec;
    }

    // ***********************************************************************//
    //                 Update Single Token Topic Assignment                   //
    // ***********************************************************************//

    arma::vec update_single_token_topic_assignment(
            arma::cube edge_probabilities,
            int tokens_in_document,
            int current_token_topic_assignment,
            arma::vec current_document_topic_counts,
            arma::mat word_type_topic_counts,
            arma::vec topic_token_counts,
            int current_word_type,
            arma::vec alpha_m,
            arma::vec beta_n,
            arma::vec document_edge_values,
            arma::vec topic_interaction_patterns,
            int document_sender,
            double rand_num,
            bool parallel,
            bool use_cached_token_topic_distribution,
            arma::vec cached_token_topic_distribution) {

        // calculate beta, the sum over beta_n
        double beta = arma::sum(beta_n);

        // get the number of topics
        int number_of_topics = alpha_m.n_elem;

        // generate a zero vector of length (number of topics)
        arma::vec unnormalized_log_distribution = arma::zeros(number_of_topics);

        // we can form the token-topic distribution in parallel since each entry
        // is independent of the others. I have verified that we get exactly the
        // same thing out either way.
        if (parallel) {
            // if we are using the cached values, then we do not need to
            // reconstruct the distribution.
            use_cached_token_topic_distribution = false;
            if (use_cached_token_topic_distribution) {
                unnormalized_log_distribution = cached_token_topic_distribution;
            } else {
                // use RcppParallel::parallelFor to loop over the distribution over
                // topics.
                unnormalized_log_distribution = parallel_token_topic_probabilities(
                    edge_probabilities,
                    tokens_in_document,
                    current_token_topic_assignment,
                    current_document_topic_counts,
                    word_type_topic_counts,
                    topic_token_counts,
                    current_word_type,
                    alpha_m,
                    beta_n,
                    document_edge_values,
                    topic_interaction_patterns,
                    document_sender,
                    beta,
                    number_of_topics);
            }

        // if we do not want to do this in parallel (which may actually be
        // faster for small numbers of topics, then just loop over the topics).
        } else {
            // loop over the topics to populate the unnormailized log distribution.
            for (int i = 0; i < number_of_topics; ++i) {
                double lsm_contr = mjd::lsm_contribution (
                    edge_probabilities,
                    tokens_in_document,
                    i,
                    current_token_topic_assignment,
                    current_document_topic_counts,
                    document_edge_values,
                    topic_interaction_patterns,
                    document_sender);
                // for testing.
                // Rcpp::Rcout << "LSM Contribution: " << lsm_contr << std::endl;

                double lda_contr = mjd::lda_contribution(
                    tokens_in_document,
                    current_token_topic_assignment,
                    current_document_topic_counts,
                    word_type_topic_counts,
                    topic_token_counts,
                    i,
                    current_word_type,
                    alpha_m,
                    beta_n,
                    beta);
                // for testing.
                // Rcpp::Rcout << "LDA Contribution: " << lda_contr << std::endl;

                // add the values (since we are working in log space) and put them
                // in the appropriate bin the in the distribution.
                unnormalized_log_distribution[i] = lsm_contr + lda_contr;
            }
        }

        // now we need to sample from this distribution
        int new_assignment = mjd::log_space_multinomial_sampler (
            unnormalized_log_distribution,
            rand_num);

        // put the token topic assignment in the first spot and the dist
        // in the rest.
        arma::vec ret_vec = arma::zeros(number_of_topics + 1);
        ret_vec[0] = new_assignment;
        for (int i = 0; i < number_of_topics; ++i) {
            ret_vec[i+1] = unnormalized_log_distribution[i];
        }

        return ret_vec;
    }


    // ***********************************************************************//
    //                 Update All Token Topic Assignments                     //
    // ***********************************************************************//

    Rcpp::List update_token_topic_assignments(
            arma::vec author_indexes,
            arma::mat document_edge_matrix,
            arma::vec topic_interaction_patterns,
            arma::mat document_topic_counts,
            arma::mat word_type_topic_counts,
            arma::vec topic_token_counts,
            Rcpp::List token_topic_assignments,
            Rcpp::List token_word_types,
            arma::vec intercepts,
            arma::mat coefficients,
            arma::cube latent_positions,
            arma::cube covariates,
            arma::vec alpha_m,
            arma::vec beta_n,
            arma::vec random_numbers,
            bool using_coefficients,
            bool parallel) {

        // get important constants
        int number_of_documents = document_edge_matrix.n_rows;
        int number_of_actors = document_edge_matrix.n_cols;
        int number_of_interaction_patterns = intercepts.n_elem;
        int rand_num_counter = 0;
        int number_of_topics = alpha_m.n_elem;

        arma::cube edge_probabilities = arma::zeros(number_of_actors,
            number_of_actors,
            number_of_interaction_patterns);

        // loop over docs, actors, interaction patterns to fill edge_probabilities
        for (int i = 0; i < number_of_actors; ++i) {
            for (int j = 0; j < number_of_actors; ++j) {
                // .tube gives us all slices
                arma::vec current_covariates = covariates.tube(i,j);
                if (i != j) {
                    for (int k = 0; k < number_of_interaction_patterns; ++k) {
                        edge_probabilities(i,j,k) = mjd::edge_probability(
                            intercepts,
                            coefficients,
                            latent_positions,
                            i,
                            j,
                            current_covariates,
                            k,
                            using_coefficients);
                    }
                }
            }
        }

        // loop over documents
        for (int i = 0; i < number_of_documents; ++i) {
            // values to be used for caching if parallel == TRUE
            bool tta_did_not_change = false;
            arma::vec cached_token_topic_distribution = arma::zeros(number_of_topics);
            bool previous_topic_assign_same_as_current = false;
            int previous_topic_assign = -1;
            bool use_cached_token_topic_distribution = false;

            // get the current token topic assignments as a vector
            arma::vec current_token_topic_assignments = token_topic_assignments[i];
            // get the current token word types as a vector
            arma::vec current_token_word_types = token_word_types[i];
            // get the current number of tokens
            int tokens_in_document = current_token_topic_assignments.n_elem;
            // allocate all of our document specific variables

            // we need to allocate row vectors from a matrix before we can
            // convert to arma::vec as desired (typeing issue)
            arma::rowvec temp = document_topic_counts.row(i);
            arma::rowvec temp2 = document_edge_matrix.row(i);

            // now we convert to an arma::vec (ugly but it works)
            arma::vec current_document_topic_counts = arma::conv_to<arma::vec>::from(temp);
            arma::vec document_edge_values = arma::conv_to<arma::vec>::from(temp2);

            int document_sender = author_indexes[i];

            // loop over tokens
            for (int j = 0; j < tokens_in_document; ++j) {
                // allocate all of our token specific variables
                int current_token_topic_assignment = current_token_topic_assignments[j];
                int current_word_type = current_token_word_types[j];
                double rand_num = random_numbers[rand_num_counter];
                rand_num_counter += 1;

                if (previous_topic_assign == current_token_topic_assignment) {
                    previous_topic_assign_same_as_current = true;
                } else {
                    previous_topic_assign_same_as_current = false;
                }

                // now determine whether we should use caching
                if (previous_topic_assign_same_as_current &
                    tta_did_not_change) {
                    use_cached_token_topic_distribution = true;
                } else {
                    use_cached_token_topic_distribution = false;
                }


                // now get the new assignment
                arma::vec new_topic_assignment_vec = update_single_token_topic_assignment(
                    edge_probabilities,
                    tokens_in_document,
                    current_token_topic_assignment,
                    current_document_topic_counts,
                    word_type_topic_counts,
                    topic_token_counts,
                    current_word_type,
                    alpha_m,
                    beta_n,
                    document_edge_values,
                    topic_interaction_patterns,
                    document_sender,
                    rand_num,
                    parallel,
                    use_cached_token_topic_distribution,
                    cached_token_topic_distribution);

                // extract values from returned vec
                int new_topic_assignment = new_topic_assignment_vec[0];
                for (int i = 0; i < number_of_topics; ++i) {
                    cached_token_topic_distribution[i] =
                        new_topic_assignment_vec[i+1];;
                }

                // set cacheing values
                tta_did_not_change = true;
                previous_topic_assign = new_topic_assignment;
                // if the assignment changed, then we need to update everything.
                if (new_topic_assignment != current_token_topic_assignment) {
                    document_topic_counts(i,current_token_topic_assignment) -= 1;
                    document_topic_counts(i,new_topic_assignment) += 1;
                    current_token_topic_assignments[j] = new_topic_assignment;
                    topic_token_counts[current_token_topic_assignment] -= 1;
                    topic_token_counts[new_topic_assignment] += 1;
                    word_type_topic_counts(current_word_type,
                                           current_token_topic_assignment) -= 1;
                    word_type_topic_counts(current_word_type,
                                           new_topic_assignment) += 1;
                    tta_did_not_change = false;
                }
            }
            // put the vector back in the list
            token_topic_assignments[i] = current_token_topic_assignments;
        }

        // put everything in a list and return it
        Rcpp::List to_return(5);
        to_return[0] = document_topic_counts;
        to_return[1] = word_type_topic_counts;
        to_return[2] = topic_token_counts;
        to_return[3] = token_topic_assignments;
        to_return[4] = edge_probabilities;

        return to_return;
    }


    // ***********************************************************************//
    //                 Update Interaction Pattern Parameters                  //
    // ***********************************************************************//

    Rcpp::List update_interaction_pattern_parameters(
            arma::vec author_indexes,
            arma::mat document_edge_matrix,
            arma::mat document_topic_counts,
            arma::vec topic_interaction_patterns,
            arma::vec intercepts,
            arma::mat coefficients,
            arma::cube latent_positions,
            arma::cube covariates,
            bool using_coefficients,
            double intercept_prior_mean,
            double intercept_prior_standard_deviation,
            arma::vec intercept_proposal_standard_deviations,
            double coefficient_prior_mean,
            double coefficient_prior_standard_deviation,
            arma::vec coefficient_proposal_standard_deviations,
            double latent_position_prior_mean,
            double latent_position_prior_standard_deviation,
            arma::vec latent_position_proposal_standard_deviations,
            double random_number,
            arma::cube edge_probabilities) {

        // get important constants
        int number_of_documents = document_edge_matrix.n_rows;
        int number_of_actors = document_edge_matrix.n_cols;
        int number_of_interaction_patterns = intercepts.n_elem;

        // draw proposed interaction parameters
        Rcpp::List proposed_parameters = mjd::sample_new_interaction_pattern_parameters(
            intercepts,
            coefficients,
            latent_positions,
            intercept_proposal_standard_deviations,
            coefficient_proposal_standard_deviations,
            latent_position_proposal_standard_deviations,
            using_coefficients);

        // allocate the appropriate objects out of the list returned above
        arma::vec proposed_intercepts = proposed_parameters[0];
        arma::mat proposed_coefficients = proposed_parameters[1];
        arma::cube proposed_latent_positions = proposed_parameters[2];

        // get the log prior probability of the current interaction pattern
        // parameters
        double log_current_prior = mjd::prior_probability_interaction_pattern_parameters(
            intercepts,
            coefficients,
            latent_positions,
            intercept_prior_mean,
            intercept_prior_standard_deviation,
            coefficient_prior_mean,
            coefficient_prior_standard_deviation,
            latent_position_prior_mean,
            latent_position_prior_standard_deviation,
            using_coefficients);

        // get the log prior probability of the proposed interaction pattern
        // parameters.
        double log_proposed_prior = mjd::prior_probability_interaction_pattern_parameters(
            proposed_intercepts,
            proposed_coefficients,
            proposed_latent_positions,
            intercept_prior_mean,
            intercept_prior_standard_deviation,
            coefficient_prior_mean,
            coefficient_prior_standard_deviation,
            latent_position_prior_mean,
            latent_position_prior_standard_deviation,
            using_coefficients);

        // allocate variables outside of loop.
        double log_current_probability = 0;
        double log_proposed_probability = 0;

        arma::cube proposed_edge_probabilities = arma::zeros(number_of_actors,
            number_of_actors,
            number_of_interaction_patterns);

        // loop over docs, actors, interaction patterns to fill edge_probabilities
        for (int i = 0; i < number_of_actors; ++i) {
            for (int j = 0; j < number_of_actors; ++j) {
                // .tube gives us all slices
                arma::vec current_covariates = covariates.tube(i,j);
                if (i != j) {
                    for (int k = 0; k < number_of_interaction_patterns; ++k) {
                        proposed_edge_probabilities(i,j,k) = mjd::edge_probability(
                            proposed_intercepts,
                            proposed_coefficients,
                            proposed_latent_positions,
                            i,
                            j,
                            current_covariates,
                            k,
                            using_coefficients);
                    }
                }
            }
        }

        // loop over documents
        for (int i = 0; i < number_of_documents; ++i) {

            // allocate all of our document specific variables

            // we need to allocate row vectors from a matrix before we can
            // convert to arma::vec as desired (typeing issue)
            arma::rowvec temp = document_topic_counts.row(i);
            arma::rowvec temp2 = document_edge_matrix.row(i);

            // now we convert to an arma::vec (ugly but it works)
            arma::vec current_document_topic_counts = arma::conv_to<arma::vec>::from(temp);
            arma::vec document_edge_values = arma::conv_to<arma::vec>::from(temp2);
            // get the current number of tokens
            int tokens_in_document = arma::sum(current_document_topic_counts);
            int document_sender = author_indexes[i];
            // loop over tokens
            // for testing
            // Rcpp::Rcout << "Document " << i << std::endl;
            for (int j = 0; j < number_of_actors; ++j) {
                // if the assignment changed, then we need to update everything.
                if (document_sender != j) {
		  double temp = sum_over_t_edge_probability(
		      edge_probabilities,
		      tokens_in_document,
		      -1,
		      current_document_topic_counts,
		      false,
		      topic_interaction_patterns,
		      document_sender,
		      j,
		     -1);

		  double temp2 = sum_over_t_edge_probability (
		      proposed_edge_probabilities,
		      tokens_in_document,
		      -1,
		      current_document_topic_counts,
		      false,
		      topic_interaction_patterns,
		      document_sender,
		      j,
		      -1);

		  if (document_edge_values[j] == 1) {
		    log_current_probability += log(temp);
		    log_proposed_probability += log(temp2);
		  } else {
		    log_current_probability += log(1 - temp);
		    log_proposed_probability += log(1 - temp2);
		  }
                } // end of condition making sure actor is not author
            } // end of loop over actors
        }// end of loop over documents

        double accept_log_prob = log_proposed_probability + log_proposed_prior -
             log_current_probability - log_current_prior;
        double log_random_number = log(random_number);
        Rcpp::List to_return(11);

        if (accept_log_prob > log_random_number) {
            to_return[0] = proposed_intercepts;
            to_return[1] = proposed_coefficients;
            to_return[2] = proposed_latent_positions;
            to_return[3] = proposed_edge_probabilities;
            to_return[4] = 1; //tells us whether we accepted
            to_return[5] = log_current_probability;
	        to_return[6] = intercepts;
	        to_return[7] = coefficients;
	        to_return[8] = latent_positions;
	        to_return[9] = edge_probabilities;
	        to_return[10] = log_proposed_probability;
        } else {
            to_return[0] = intercepts;
            to_return[1] = coefficients;
            to_return[2] = latent_positions;
            to_return[3] = edge_probabilities;
            to_return[4] = 0; //tells us whether we accepted
            to_return[5] = log_current_probability;
	        to_return[6] = proposed_intercepts;
            to_return[7] = proposed_coefficients;
            to_return[8] = proposed_latent_positions;
            to_return[9] = proposed_edge_probabilities;
	        to_return[10] = log_proposed_probability;
        }

        return to_return;
    }


    // ***********************************************************************//
    //                 Update Topic Interaction Pattern Assignments           //
    // ***********************************************************************//

    arma::vec update_topic_interaction_pattern_assignments(
            arma::vec author_indexes,
            arma::mat document_edge_matrix,
            arma::mat document_topic_counts,
            arma::vec topic_interaction_patterns,
            arma::vec intercepts,
            arma::mat coefficients,
            arma::cube latent_positions,
            arma::cube covariates,
            bool using_coefficients,
            arma::vec random_numbers,
            arma::cube edge_probabilities) {

        // get important constants
        int number_of_documents = document_edge_matrix.n_rows;
        int number_of_actors = document_edge_matrix.n_cols;
        int number_of_interaction_patterns = intercepts.n_elem;
        int number_of_topics = topic_interaction_patterns.n_elem;
        int number_counter = 0;

        // outer loop over topics
        for (int t = 0; t < number_of_topics; ++t) {
            //allocate topic specific variables
            arma::vec interaction_pattern_assignment_log_probs = arma::zeros(
                number_of_interaction_patterns);
            arma::mat held_out_sum_over_t_terms = arma::zeros(
                number_of_documents, number_of_actors);

            //calculate held out sum over t terms
            // loop over documents
            for (int i = 0; i < number_of_documents; ++i) {
                // allocate all of our document specific variables

                // we need to allocate row vectors from a matrix before we can
                // convert to arma::vec as desired (typeing issue)
                arma::rowvec temp = document_topic_counts.row(i);
                arma::rowvec temp2 = document_edge_matrix.row(i);

                // now we convert to an arma::vec (ugly but it works)
                arma::vec current_document_topic_counts = arma::conv_to<arma::vec>::from(temp);
                arma::vec document_edge_values = arma::conv_to<arma::vec>::from(temp2);

                // get the current number of tokens
                int tokens_in_document = arma::sum(current_document_topic_counts);
                int document_sender = author_indexes[i];
                // loop over tokens
                for (int j = 0; j < number_of_actors; ++j) {
                    // if the assignment changed, then we need to update everything.
                    if (document_sender != j) {
		      double temp = sum_over_t_edge_probability (
                          edge_probabilities,
                          tokens_in_document,
                          -1,
                          current_document_topic_counts,
                          false,
                          topic_interaction_patterns,
                          document_sender,
                          j,
                          t);

                        if (document_edge_values[j] == 1) {
                            held_out_sum_over_t_terms(i,j) = temp;
                        } else {
                            held_out_sum_over_t_terms(i,j) = 1 - temp;
                        }
                    } // end of condition making sure actor is not author
                } // end of loop over actors
            }// end of loop over documents

            // loop over interaction patterns to populate distribution
            for (int c = 0; c < number_of_interaction_patterns; ++c) {
                double log_prob = 0;
                // loop over documents
                for (int i = 0; i < number_of_documents; ++i) {
                    // allocate all of our document specific variables

                    // we need to allocate row vectors from a matrix before we can
                    // convert to arma::vec as desired (typeing issue)
                    arma::rowvec temp = document_topic_counts.row(i);
                    arma::rowvec temp2 = document_edge_matrix.row(i);

                    // now we convert to an arma::vec (ugly but it works)
                    arma::vec current_document_topic_counts = arma::conv_to<arma::vec>::from(temp);
                    arma::vec document_edge_values = arma::conv_to<arma::vec>::from(temp2);

                    // get the current number of tokens
                    int tokens_in_document = arma::sum(current_document_topic_counts);
                    int document_sender = author_indexes[i];
                    // loop over tokens
                    for (int j = 0; j < number_of_actors; ++j) {
                        // if the assignment changed, then we need to update everything.
                        if (document_sender != j) {
                            double ttc = double(current_document_topic_counts[t])/
                            double(tokens_in_document);
                            if (document_edge_values[j] == 1) {
                                log_prob += log(held_out_sum_over_t_terms(i,j) +
                                    double(ttc *
                                    edge_probabilities(document_sender,j,c)));
                            } else {
                                log_prob += log(held_out_sum_over_t_terms(i,j) +
                                    1 - double(ttc *
                                    edge_probabilities(document_sender,j,c)));
                            }
                        } // end of condition making sure actor is not author
                    } // end of loop over actors
                }// end of loop over documents
                interaction_pattern_assignment_log_probs[c] = log_prob;
            }//end of loop over clusters
            double random_number = random_numbers[number_counter];
            number_counter += 1;
            int new_assignment = log_space_multinomial_sampler(
                interaction_pattern_assignment_log_probs,
                random_number);

            topic_interaction_patterns[t] =  new_assignment;
        }//end of loop over topics

        //return
        return topic_interaction_patterns;
    }

    // ***********************************************************************//
    //                             Adaptive Metropolis                        //
    // ***********************************************************************//

    Rcpp::List adaptive_metropolis(
            arma::vec intercept_proposal_standard_deviations,
            arma::vec coefficient_proposal_standard_deviations,
            arma::vec latent_position_proposal_standard_deviations,
            arma::vec accept_rates,
            double target_accept_rate,
            double tollerance,
            double update_size) {

        Rcpp::Rcout << "Metropolis Hastings Accept Rates: " << std::endl <<
            accept_rates.t() << std::endl;

        // get number of interaction patterns
        int number_of_interaction_patterns = intercept_proposal_standard_deviations.n_elem;

        // loop over interaction patterns (in the future, we could have different
        // accept rates for different interaction patterns).
        for (int i = 0; i < number_of_interaction_patterns; ++i) {
            // if our accept rate is too low, bump down proposal sd
            if (accept_rates[i] < (target_accept_rate - tollerance)) {
                // for now, we just want to make sure we do not go below zero
                double floor = intercept_proposal_standard_deviations[i] - update_size;

                // only update if it will not make proposal standard_deviation negative
                if (floor > 0) {
                    intercept_proposal_standard_deviations[i] -= update_size;
                    coefficient_proposal_standard_deviations[i] -= update_size;
                    latent_position_proposal_standard_deviations[i] -= update_size;
                }
            }
            // if the accept rate is too high, bump up the proposal sd
            if (accept_rates[i] > (target_accept_rate + tollerance)) {
                intercept_proposal_standard_deviations[i] += update_size;
                coefficient_proposal_standard_deviations[i] += update_size;
                latent_position_proposal_standard_deviations[i] += update_size;
            }
        }

        // print out the new proposal sd's so the user can keep track of what
        // is going on.
        Rcpp::Rcout << "New Proposal Standard Deviations: " << std::endl <<
            intercept_proposal_standard_deviations.t() << std::endl;

        // put everything in a return list
        Rcpp::List ret_list(3);
        ret_list[0] = intercept_proposal_standard_deviations;
        ret_list[1] = coefficient_proposal_standard_deviations;
        ret_list[2] = latent_position_proposal_standard_deviations;

        return ret_list;
    }

    // ***********************************************************************//
    //           Calculate the unnormalized LDA log likelihood                //
    // ***********************************************************************//

    // this code is adapted from from Hanna's Github -- need to get exact source
    double calculate_unnormalized_LDA_log_likelihood(int number_of_documents,
                                                     int number_of_topics,
                                                     Rcpp::List token_topic_assignments,
                                                     arma::vec log_alpha_m) {

        //initialize variables that will be used across iterations
        double Unnormalized_Corpus_Log_Likelihood = 0;

        for(int d = 0; d < number_of_documents; ++d){
            //set all document specific parameters
            arma::vec current_token_topic_assignments = token_topic_assignments[d];
            int number_of_tokens = current_token_topic_assignments.n_elem;
            double current_document_contribution = 0;
            double current_document_contribution2 = 0;

            //loop over topics
            for(int t = 0; t < number_of_topics; ++t){
                int topic = t + 1;
                // now we calculate the first and second terms in the likelihood of
                // the token being from the current topic calculate the number of
                // times a token in the current document has been assigned to the
                // current topic.
                int ntd = 0;
                for(int b = 0; b < number_of_tokens; ++b){
                    if(current_token_topic_assignments[b] == topic){
                        ntd +=1;
                    }
                }
                double first_term = double(ntd) + exp(log_alpha_m[t]);
                current_document_contribution += lgamma(first_term);
                double first_term2 = exp(log_alpha_m[t]);
                current_document_contribution2 -= lgamma(first_term2);
            }//loop over topics

            double second_term = number_of_tokens +
                exp(double(number_of_topics)*log_alpha_m[0]);
            current_document_contribution -= lgamma(second_term);

            double second_term2 = exp(double(number_of_topics)*log_alpha_m[0]);
            current_document_contribution2 += second_term2;

            Unnormalized_Corpus_Log_Likelihood += current_document_contribution +
                current_document_contribution2;

        }//end of document loop

        // we have to add on an alpha at the end because we are doing the log
        // transform
        Unnormalized_Corpus_Log_Likelihood += (double(number_of_topics) *
            log_alpha_m[0]);

        Rcpp::Rcout << "Corpus Log Likelihood: " <<
            Unnormalized_Corpus_Log_Likelihood << std::endl;

        return Unnormalized_Corpus_Log_Likelihood;
    }

    // ***********************************************************************//
    //                  Slice sample new values for alpha_m                   //
    // ***********************************************************************//

    arma::vec slice_sample_new_alpha_m(int number_of_documents,
                                    int number_of_topics,
                                    Rcpp::List token_topic_assignments,
                                    arma::vec log_alpha_m,
                                    double current_LDA_ll,
                                    double step_size) {

        // use R::runif becasue we are using a while loop so we do not know how many
        // numbers we will need to pass in. It is slower but safer.
        double rand_num = R::runif(0,1);
        double lud = log(rand_num);
        double slice_probability_floor = current_LDA_ll + lud;

        arma::vec proposed_alpha_m = arma::zeros(number_of_topics);
        arma::vec left_proposed_alpha_m = arma::zeros(number_of_topics);
        arma::vec right_proposed_alpha_m = arma::zeros(number_of_topics);
        //get the left and right bounds on the slice
        rand_num = R::runif(0,1);
        for(int t = 0; t < number_of_topics; ++t){
            left_proposed_alpha_m[t] = log_alpha_m[t] - rand_num*step_size;
            right_proposed_alpha_m[t] = left_proposed_alpha_m[t] + step_size;
        }

        //set equal to one when new sample accepted
        int in_slice = 0;

        while(in_slice < 1){
            //report(proposed_alpha_m[0]);
            //get new values for the slice for alpha
            rand_num = R::runif(0,1);
            for(int t = 0; t < number_of_topics; ++t){
                proposed_alpha_m[t] = left_proposed_alpha_m[t] +
                    rand_num*(right_proposed_alpha_m[t] - left_proposed_alpha_m[t]);
            }

            double proposed_LDA_ll = calculate_unnormalized_LDA_log_likelihood(
                number_of_documents,
                number_of_topics,
                token_topic_assignments,
                proposed_alpha_m);

            // ========== check to see if it is under the curve ======== //
            if(proposed_LDA_ll > slice_probability_floor){
                in_slice = 1;
            }
            else{
                //if the positions we tried were outside of the slice, set them as the new boundary
                //get the left and right bounds on the slice for alpha
                for(int t = 0; t < number_of_topics; ++t){
                    if(proposed_alpha_m[t] < log_alpha_m[t]){
                        left_proposed_alpha_m[t] = proposed_alpha_m[t];
                    }
                    else{
                        right_proposed_alpha_m[t] = proposed_alpha_m[t];
                    }
                }
            }
        }// end of while checking to see if we are in slice loop

        arma::vec alpha_m = arma::zeros(number_of_topics);
        for(int t = 0; t < number_of_topics; ++t){
            alpha_m[t] = exp(proposed_alpha_m[t]);
        }
        return alpha_m;
    }

    // ***********************************************************************//
    //                             Inference                                  //
    // ***********************************************************************//

    Rcpp::List inference(
            arma::vec author_indexes,
            arma::mat document_edge_matrix,
            arma::mat document_topic_counts,
            arma::vec topic_interaction_patterns,
            arma::mat word_type_topic_counts,
            arma::vec topic_token_counts,
            Rcpp::List token_topic_assignments,
            Rcpp::List token_word_types,
            arma::vec intercepts,
            arma::mat coefficients,
            arma::cube latent_positions,
            arma::cube covariates,
            arma::vec alpha_m,
            arma::vec beta_n,
            bool using_coefficients,
            double intercept_prior_mean,
            double intercept_prior_standard_deviation,
            arma::vec intercept_proposal_standard_deviations,
            double coefficient_prior_mean,
            double coefficient_prior_standard_deviation,
            arma::vec coefficient_proposal_standard_deviations,
            double latent_position_prior_mean,
            double latent_position_prior_standard_deviation,
            arma::vec latent_position_proposal_standard_deviations,
            double target_accept_rate,
            double tollerance,
            double update_size,
            int seed,
            int iterations,
            int metropolis_iterations,
            int total_number_of_tokens,
            int iterations_before_t_i_p_updates,
            int update_t_i_p_every_x_iterations,
            bool perform_adaptive_metropolis,
            int slice_sample_every_x_iterations,
            double slice_sample_step_size,
            bool parallel) {

        // Set RNG and define uniform distribution
        boost::mt19937 generator(seed);
        boost::uniform_01<double> uniform_distribution;

        // example get the random uniform draw
        // double rand_num = uniform_distribution(generator);

        // allocated global variables
        int t_i_p_update_counter = 0;
        int slice_sample_counter = 0;
        int num_interaction_patterns = intercept_proposal_standard_deviations.n_elem;
        int num_topics = topic_interaction_patterns.n_elem;
        int num_actors = document_edge_matrix.n_cols;
        int num_latent_dimensions = latent_positions.n_slices;
        int number_of_documents = document_edge_matrix.n_rows;
        //get the number of if we are using them, set this eqaul to two since we
        //need to allocate the matrix even if we are not using it.
        int num_coefficients = 2;
        if(using_coefficients){
            num_coefficients = coefficients.n_cols;
        }

        arma::cube edge_probabilities = arma::zeros(num_actors, num_actors,
                                              num_interaction_patterns);
        arma::vec accept_rates = arma::zeros(num_interaction_patterns);

        // allocate data structures to store samples in.
        arma::mat store_topic_interaction_patterns = arma::zeros(
            iterations, num_topics);
        arma::mat store_intercepts = arma::zeros(metropolis_iterations,
            num_interaction_patterns);
        arma::cube store_coefficients = arma::zeros(num_interaction_patterns,
            num_coefficients, metropolis_iterations);
        //we are going to have to stack cubes here, hence the multiplication
        //in the last dimension
        arma::cube store_latent_positions = arma::zeros(
            num_interaction_patterns,
            num_actors,
            metropolis_iterations * num_latent_dimensions);

        arma::mat store_accept_rates = arma::zeros(iterations,
                                                   num_interaction_patterns);
        arma::vec store_LDA_ll = arma::zeros(iterations);

        // loop over interaction patterns
        for (int i = 0; i < iterations; ++i) {
            Rcpp::Rcout << "Iteration: " << i << std::endl;

            // generate a vector of random numbers to pass in to the topic-token
            // update function.
            arma::vec random_numbers = arma::zeros(total_number_of_tokens);
            for (int k = 0; k < total_number_of_tokens; ++k) {
                random_numbers[k] = uniform_distribution(generator);
            }

            // update all token topic assignments
            Rcpp::List Topic_Updates = update_token_topic_assignments(
                author_indexes,
                document_edge_matrix,
                topic_interaction_patterns,
                document_topic_counts,
                word_type_topic_counts,
                topic_token_counts,
                token_topic_assignments,
                token_word_types,
                intercepts,
                coefficients,
                latent_positions,
                covariates,
                alpha_m,
                beta_n,
                random_numbers,
                using_coefficients,
                parallel);

            // take everything out of the returned list and put it back in the
            // main objects. We have to do this silly double assignment because
            // the direct assignment from and Rcpp::List object is ambiguous.
            arma::mat temp = Topic_Updates[0];
            document_topic_counts = temp;
            arma::mat temp2 = Topic_Updates[1];
            word_type_topic_counts = temp2;
            arma::vec temp3 = Topic_Updates[2];
            topic_token_counts = temp3;
            //can do direct assignment here because the list entry is itself a
            //list object.
            token_topic_assignments = Topic_Updates[3];
            arma::cube temp4 = Topic_Updates[4];
            edge_probabilities = temp4;

            // only update topic interaction pattern assignments if we have
            // completed atleast x iterations.
            if (i >= iterations_before_t_i_p_updates) {
                // this conditional allows us to only update topic interaction
                // pattern assignments every x iterations
                if (t_i_p_update_counter >= update_t_i_p_every_x_iterations) {

                    // generate a vector of random numbers to pass in to the topic-token
                    // update function.
                    arma::vec random_numbers = arma::zeros(num_topics);
                    for (int k = 0; k < num_topics; ++k) {
                        random_numbers[k] = uniform_distribution(generator);
                    }
                    //update topic interaction patterns
                    topic_interaction_patterns = update_topic_interaction_pattern_assignments(
                        author_indexes,
                        document_edge_matrix,
                        document_topic_counts,
                        topic_interaction_patterns,
                        intercepts,
                        coefficients,
                        latent_positions,
                        covariates,
                        using_coefficients,
                        random_numbers,
                        edge_probabilities);
                    //reset counter
                    t_i_p_update_counter = 0;

                    // print the T. I. P. Assignments out so we can check for weird
                    // stuff like all topics going toone IP.
                    Rcpp::Rcout << "Topic Interaction Pattern Assignments: " <<
                        std::endl << Rcpp::as<Rcpp::NumericVector>(
                        Rcpp::wrap(topic_interaction_patterns.t())) << std::endl;
                }
                //increment counter
                t_i_p_update_counter += 1;
            } // end of topic interaction pattern update loop

            // ************************************************************** //
            // This section of the code is not necessary for our algorithm to
            // perform inference, but will give us convergence diagnostics and
            // better fit

            // we need to work with the log of alpha m for slice sampling
            arma::vec log_alpha_m(num_topics);
            for(int t = 0; t < num_topics; ++t){
                log_alpha_m[t] = log(alpha_m[t]);
            }

            // get the unnormalied LDA log likelihood for assesing convergence.
            double LDA_ll = calculate_unnormalized_LDA_log_likelihood(
                number_of_documents,
                num_topics,
                token_topic_assignments,
                log_alpha_m);

            // store the value.
            store_LDA_ll[i] = LDA_ll;

            // if we are slice sampling this iteration, the do it.
            if (slice_sample_every_x_iterations == slice_sample_counter) {
                alpha_m = slice_sample_new_alpha_m(number_of_documents,
                                                   num_topics,
                                                   token_topic_assignments,
                                                   log_alpha_m,
                                                   LDA_ll,
                                                   slice_sample_step_size);
                slice_sample_counter = 0;
                Rcpp::Rcout << "New Alpha value:" << alpha_m[0] << std::endl;
            }
            slice_sample_counter += 1;
            // ************************************************************** //

            // perform adaptive metropolis updates if there has been atleast
            // one outer iteration
            if (i > 0) {
                // only perform adaptive metropolis if it is turned on via the
                // boolean
                if (perform_adaptive_metropolis) {
                    Rcpp::List am_updates =  adaptive_metropolis(
                            intercept_proposal_standard_deviations,
                            coefficient_proposal_standard_deviations,
                            latent_position_proposal_standard_deviations,
                            accept_rates,
                            target_accept_rate,
                            tollerance,
                            update_size);

                    // deal with type ambiguity while extracting objects from an
                    // Rcpp List, then assign everything to the appropriate
                    // vectors.
                    arma::vec temp = am_updates[0];
                    intercept_proposal_standard_deviations = temp;
                    arma::vec temp2 = am_updates[1];
                    coefficient_proposal_standard_deviations = temp2;
                    arma::vec temp3 = am_updates[2];
                    latent_position_proposal_standard_deviations = temp3;

                }// end of conditional for whether we perform adaptive MH
            }// end of conditional checking to see if we are past the first update.

            // these variables will be used for storage
            int latent_position_storage_counter = 0;
            arma::vec accept_or_reject = arma::zeros(metropolis_iterations);
            // loop over metropolis hastings iterations
            for (int j = 0; j < metropolis_iterations; ++j) {
                double random_number =  uniform_distribution(generator);
                Rcpp::List MH_List =  update_interaction_pattern_parameters(
                        author_indexes,
                        document_edge_matrix,
                        document_topic_counts,
                        topic_interaction_patterns,
                        intercepts,
                        coefficients,
                        latent_positions,
                        covariates,
                        using_coefficients,
                        intercept_prior_mean,
                        intercept_prior_standard_deviation,
                        intercept_proposal_standard_deviations,
                        coefficient_prior_mean,
                        coefficient_prior_standard_deviation,
                        coefficient_proposal_standard_deviations,
                        latent_position_prior_mean,
                        latent_position_prior_standard_deviation,
                        latent_position_proposal_standard_deviations,
                        random_number,
                        edge_probabilities);

                // extract everything from the list returned by the MH parameter
                // update function.
                arma::vec temp = MH_List[0];
                intercepts = temp;
                arma::mat temp2 = MH_List[1];
                coefficients = temp2;
                arma::cube temp3 = MH_List[2];
                latent_positions = temp3;
                arma::cube temp4 = MH_List[3];
                edge_probabilities = temp4;
                Rcpp::NumericVector temp5a = MH_List[4];
                double temp5 = temp5a[0];
                accept_or_reject[j] = temp5;

                // if we are on the last iteration of gibbs sampling, then we
                // need to save everything
                if (i == (iterations-1)) {
                   for (int k = 0; k < num_interaction_patterns; ++k) {
                       //save intercepts
                       double temp = intercepts[k];
                       store_intercepts(j,k) = temp;

                       if (using_coefficients) {
                           for (int l = 0; l < num_coefficients; ++l) {
                               double temp2 = coefficients(k,l);
                               store_coefficients(k,l,j) = temp2;
                           }
                       }
                   }
                   // store latent positions -- we need the latent dimensions
                   // to be the outter loop so that we can increment the counter
                   // of which slice we insert into.
                   for (int m = 0; m < num_latent_dimensions; ++m) {
                       for (int k = 0; k < num_interaction_patterns; ++k) {
                           for (int l = 0; l < num_actors; ++l) {
                               double temp3 = latent_positions(k,l,m);
                               store_latent_positions(k,l,
                                    latent_position_storage_counter) = temp3;
                           }
                       }
                       latent_position_storage_counter += 1;
                   }// loop over latent dimensions has to be last as they are
                   // stacked in the returned array
                }//end of storage conditional

            }// end of metropolis hastings loop

            // store topic interaction pattern assignments
            for (int k = 0; k < num_topics; ++k) {
                int temp = topic_interaction_patterns[k];
                store_topic_interaction_patterns(i,k) = temp;
            }
            // store MCMC accept rates for each iteration
            for (int k = 0; k < num_interaction_patterns; ++k) {
                // calculate the accept proportion
                double temp = arma::sum(accept_or_reject)/double(
                    metropolis_iterations);
                // in this case, we are constraining it to be the same across all
                // interaction patterns
                accept_rates[k] = temp;
                store_accept_rates(i,k) = temp;
            }

        }// end up gibbs/main sampling loop

        // allocate a list to store everything in.
        Rcpp::List ret_list(16);
        ret_list[0] = store_topic_interaction_patterns;
        ret_list[1] = store_intercepts;
        ret_list[2] = store_coefficients;
        ret_list[3] = store_latent_positions;
        ret_list[4] = document_topic_counts;
        ret_list[5] = word_type_topic_counts;
        ret_list[6] = topic_token_counts;
        ret_list[7] = token_topic_assignments;
        ret_list[8] = store_accept_rates;
        ret_list[9] = intercept_proposal_standard_deviations;
        ret_list[10] = store_LDA_ll;
        ret_list[11] = intercepts;
        ret_list[12] = coefficients;
        ret_list[13] = latent_positions;
        ret_list[14] = topic_interaction_patterns;
        ret_list[15] = word_type_topic_counts;
        // return everything
        return ret_list;
    }

    // ***********************************************************************//
    //             Single Draw from a Dirichlet Distribution                  //
    // ***********************************************************************//

    arma::vec rdirichlet(arma::vec alpha_m) {
        // this example is drawn from:
        // https://en.wikipedia.org/wiki/Dirichlet_distribution#Random_number_generation

        int distribution_size = alpha_m.n_elem;
        arma::vec distribution = arma::zeros(distribution_size);

        // loop through the distribution and draw Gamma variables
        for (int i = 0; i < distribution_size; ++i) {
            double cur = R::rgamma(alpha_m[i],1.0);
            distribution[i] = cur;
        }

        double gamma_sum = arma::sum(distribution);

        for (int i = 0; i < distribution_size; ++i) {
            double temp = distribution[i];
            distribution[i] = temp/gamma_sum;
        }

        return distribution;
    }

    // ***********************************************************************//
    //             Multiple Draws from a Dirichlet Distribution               //
    // ***********************************************************************//

    arma::mat rdirichlet_matrix(int num_samples,
                             arma::vec alpha_m) {
        int distribution_size = alpha_m.n_elem;
        // each row will be a draw from a Dirichlet
        arma::mat distribution = arma::zeros(num_samples, distribution_size);

        for (int i = 0; i < num_samples; ++i) {
            double sum_term = 0;
            // loop through the distribution and draw Gamma variables
            for (int j = 0; j < distribution_size; ++j) {
                double cur = R::rgamma(alpha_m[j],1.0);
                distribution(i,j) = cur;
                sum_term += cur;
            }
            // now normalize
            for (int j = 0; j < distribution_size; ++j) {
                distribution(i,j) = distribution(i,j)/sum_term;
            }
        }
        return distribution;
    }

    // ***********************************************************************//
    //        Sample Token Topic Assignments From Generative Process          //
    // ***********************************************************************//

    Rcpp::List sample_token_topics_generative_process(
            Rcpp::List token_topic_assignments,
            Rcpp::List token_word_types,
            arma::vec alpha_m,
            arma::vec beta_n,
            int number_of_documents,
            bool resample_word_types,
            arma::vec random_numbers) {

        // get some global variables
        int num_topics = alpha_m.n_elem;
        int num_word_types = beta_n.n_elem;
        int rand_num_counter = 0;

        arma::mat document_topic_counts = arma::zeros(number_of_documents,
                                                      num_topics);

        arma::vec topic_token_counts = arma::zeros(num_topics);

        arma::mat word_type_topic_counts = arma::zeros(num_word_types,
                                                       num_topics);

        // first we draw topic-word type distributions regardless of whether we
        // are going to use them
        arma::mat topic_word_type_distributions = mjd::rdirichlet_matrix(
            num_topics,
            beta_n);

        // now we draw the documen-topic distributions
        arma::mat document_topic_distributions = mjd::rdirichlet_matrix(
            number_of_documents,
            alpha_m);

        for (int i = 0; i < number_of_documents; ++i) {
            // get the current token topic assignments as a vector
            arma::vec current_token_topic_assignments = token_topic_assignments[i];
            // get the current token word types as a vector
            arma::vec current_token_word_types = token_word_types[i];
            // get the current number of tokens
            int tokens_in_document = current_token_topic_assignments.n_elem;

            // we need to allocate row vectors from a matrix before we can
            // convert to arma::vec as desired (typeing issue)
            arma::rowvec temp = document_topic_distributions.row(i);
            // now we convert to an arma::vec (ugly but it works)
            arma::vec current_doc_topic_dist = arma::conv_to<arma::vec>::from(temp);
            // need to take log so we can use our log space multinomial sampler
            current_doc_topic_dist = arma::log(current_doc_topic_dist);

            // loop over tokens
            for (int j = 0; j < tokens_in_document; ++j) {
                double rand_num = random_numbers[rand_num_counter];
                rand_num_counter += 1;
                int current_word_type = current_token_word_types[j];

                // now get the new assignment
                int new_topic_assignment = mjd::log_space_multinomial_sampler(
                    current_doc_topic_dist,
                    rand_num);

                if (resample_word_types) {
                    arma::rowvec temp2 = topic_word_type_distributions.row(new_topic_assignment);
                    // now we convert to an arma::vec (ugly but it works)
                    arma::vec current_topic_word_type_dist = arma::conv_to<arma::vec>::from(temp2);
                    // need to take log so we can use our log space multinomial sampler
                    current_topic_word_type_dist = arma::log(current_topic_word_type_dist);
                    double rand_num2 = random_numbers[rand_num_counter];
                    rand_num_counter += 1;
                    current_word_type = mjd::log_space_multinomial_sampler(
                        current_topic_word_type_dist,
                        rand_num2);
                }
                document_topic_counts(i,new_topic_assignment) += 1;
                current_token_topic_assignments[j] = new_topic_assignment;
                topic_token_counts[new_topic_assignment] += 1;
                word_type_topic_counts(current_word_type,
                                       new_topic_assignment) += 1;
                current_token_word_types[j] = current_word_type;

            }
            // put the vector back in the list
            token_topic_assignments[i] = current_token_topic_assignments;
            token_word_types[i] = current_token_word_types;
        } // end of loop over documents

        //return everything
        Rcpp::List to_return(7);
        to_return[0] = token_topic_assignments;
        to_return[1] = token_word_types;
        to_return[2] = document_topic_counts;
        to_return[3] = topic_token_counts;
        to_return[4] = word_type_topic_counts;
        to_return[5] = document_topic_distributions;
        to_return[6] = topic_word_type_distributions;

        return to_return;
    }


    // ***********************************************************************//
    //                           Run MH to Convergence                        //
    // ***********************************************************************//

    Rcpp::List run_metropolis_hastings_to_convergence(
            arma::vec author_indexes,
            arma::mat document_edge_matrix,
            arma::mat document_topic_counts,
            arma::vec topic_interaction_patterns,
            arma::vec intercepts,
            arma::mat coefficients,
            arma::cube latent_positions,
            arma::cube covariates,
            bool using_coefficients,
            double intercept_prior_mean,
            double intercept_prior_standard_deviation,
            arma::vec intercept_proposal_standard_deviations,
            double coefficient_prior_mean,
            double coefficient_prior_standard_deviation,
            arma::vec coefficient_proposal_standard_deviations,
            double latent_position_prior_mean,
            double latent_position_prior_standard_deviation,
            arma::vec latent_position_proposal_standard_deviations,
            double target_accept_rate,
            double tollerance,
            double update_size,
            int seed,
            int metropolis_iterations,
            int adaptive_metropolis_every_x_iterations,
            int stop_adaptive_metropolis_after_x_updates,
            int samples_to_store,
            int sample_every,
            int burnin) {

        // Set RNG and define uniform distribution
        boost::mt19937 generator(seed);
        boost::uniform_01<double> uniform_distribution;

        // example get the random uniform draw
        // double rand_num = uniform_distribution(generator);

        // allocated global variables
        int num_interaction_patterns = intercept_proposal_standard_deviations.n_elem;
        int num_actors = document_edge_matrix.n_cols;
        int num_latent_dimensions = latent_positions.n_slices;
        int adaptive_MH_counter = 0;
        int total_adaptive_MH_updates = 0;
        int latent_position_storage_counter = 0;
        int storage_counter = 0;
        int burnin_counter = 0;
        int storage_check = 1;
        //get the number of if we are using them, set this eqaul to two since we
        //need to allocate the matrix even if we are not using it.
        int num_coefficients = 2;
        if(using_coefficients){
            num_coefficients = coefficients.n_cols;
        }
        arma::vec accept_rates = arma::zeros(num_interaction_patterns);
        arma::vec store_accept_or_reject = arma::zeros(metropolis_iterations);
        arma::vec LSM_log_likelihood = arma::zeros(samples_to_store);
        arma::vec accept_or_reject = arma::zeros(adaptive_metropolis_every_x_iterations);

        // allocate data structures to store samples in.
        arma::mat store_intercepts = arma::zeros(samples_to_store,
                                                 num_interaction_patterns);
        arma::cube store_coefficients = arma::zeros(num_interaction_patterns,
                                                    num_coefficients,
                                                    samples_to_store);
        //we are going to have to stack cubes here, hence the multiplication
        //in the last dimension
        arma::cube store_latent_positions = arma::zeros(
            num_interaction_patterns,
            num_actors,
            samples_to_store * num_latent_dimensions);

        // generate edge probabilities to pass in
        arma::cube edge_probabilities = arma::zeros(num_actors,
                                                    num_actors,
                                                    num_interaction_patterns);

        // loop over docs, actors, interaction patterns to fill edge_probabilities
        for (int i = 0; i < num_actors; ++i) {
            for (int j = 0; j < num_actors; ++j) {
                // .tube gives us all slices
                arma::vec current_covariates = covariates.tube(i,j);
                if (i != j) {
                    for (int k = 0; k < num_interaction_patterns; ++k) {
                        edge_probabilities(i,j,k) = mjd::edge_probability(
                            intercepts,
                            coefficients,
                            latent_positions,
                            i,
                            j,
                            current_covariates,
                            k,
                            using_coefficients);
                    }
                }
            }
        }

        // loop over metropolis hastings iterations
        for (int j = 0; j < metropolis_iterations; ++j) {

            double random_number =  uniform_distribution(generator);
            Rcpp::List MH_List =  update_interaction_pattern_parameters(
                author_indexes,
                document_edge_matrix,
                document_topic_counts,
                topic_interaction_patterns,
                intercepts,
                coefficients,
                latent_positions,
                covariates,
                using_coefficients,
                intercept_prior_mean,
                intercept_prior_standard_deviation,
                intercept_proposal_standard_deviations,
                coefficient_prior_mean,
                coefficient_prior_standard_deviation,
                coefficient_proposal_standard_deviations,
                latent_position_prior_mean,
                latent_position_prior_standard_deviation,
                latent_position_proposal_standard_deviations,
                random_number,
                edge_probabilities);

            // extract everything from the list returned by the MH parameter
            // update function.
            arma::vec temp = MH_List[0];
            intercepts = temp;
            arma::mat temp2 = MH_List[1];
            coefficients = temp2;
            arma::cube temp3 = MH_List[2];
            latent_positions = temp3;
            arma::cube temp4 = MH_List[3];
            edge_probabilities = temp4;
            Rcpp::NumericVector temp5a = MH_List[4];
            double temp5 = temp5a[0];
            accept_or_reject[adaptive_MH_counter] = temp5;
            store_accept_or_reject[j] = temp5;
            Rcpp::NumericVector temp6a = MH_List[5];
            double temp6 = temp6a[0];

            // do adaptive metropolis
            if (adaptive_MH_counter >= adaptive_metropolis_every_x_iterations) {
                // print out the current iteration every x updates so we can
                // keep track of progress.
                Rcpp::Rcout << "Iteration: " << j << std::endl;
                if (total_adaptive_MH_updates <
                    stop_adaptive_metropolis_after_x_updates) {

                    for (int k = 0; k < num_interaction_patterns; ++k) {
                        // calculate the accept proportion
                        double temp = arma::sum(accept_or_reject)/double(
                            adaptive_metropolis_every_x_iterations);
                        // in this case, we are constraining it to be the same across all
                        // interaction patterns
                        accept_rates[k] = temp;
                    }

                    Rcpp::List am_updates =  adaptive_metropolis(
                        intercept_proposal_standard_deviations,
                        coefficient_proposal_standard_deviations,
                        latent_position_proposal_standard_deviations,
                        accept_rates,
                        target_accept_rate,
                        tollerance,
                        update_size);

                    // deal with type ambiguity while extracting objects from an
                    // Rcpp List, then assign everything to the appropriate
                    // vectors.
                    arma::vec temp = am_updates[0];
                    intercept_proposal_standard_deviations = temp;
                    arma::vec temp2 = am_updates[1];
                    coefficient_proposal_standard_deviations = temp2;
                    arma::vec temp3 = am_updates[2];
                    latent_position_proposal_standard_deviations = temp3;

                    // increment the total number of updates we do before
                    // stopping
                    total_adaptive_MH_updates += 1;

                }// end of conditional for whether we perform adaptive MH
                adaptive_MH_counter = 0;
            }

            // only save stuff if we are past the burnin
            if (burnin_counter > burnin) {
                // only save once every sample_every iterations
                if (storage_check == sample_every) {
                    // save everything
                    LSM_log_likelihood[storage_counter] = temp6;
                    for (int k = 0; k < num_interaction_patterns; ++k) {
                        //save intercepts
                        double temp = intercepts[k];
                        store_intercepts(storage_counter,k) = temp;

                        if (using_coefficients) {
                            for (int l = 0; l < num_coefficients; ++l) {
                                double temp2 = coefficients(k,l);
                                store_coefficients(k,l,storage_counter) = temp2;
                            }
                        }
                    }
                    // store latent positions -- we need the latent dimensions
                    // to be the outter loop so that we can increment the counter
                    // of which slice we insert into.
                    for (int m = 0; m < num_latent_dimensions; ++m) {
                        for (int k = 0; k < num_interaction_patterns; ++k) {
                            for (int l = 0; l < num_actors; ++l) {
                                double temp3 = latent_positions(k,l,m);
                                store_latent_positions(k,l,
                                                       latent_position_storage_counter) = temp3;
                            }
                        }
                        latent_position_storage_counter += 1;
                    }// loop over latent dimensions has to be last as they are
                    // stacked in the returned array

                    // this is the counter for where we stick stuff in
                    storage_counter += 1;
                    storage_check = 0;
                }

                // only increment the storage check if we have hit the burnin
                storage_check += 1;
            }// burnin conditional

            // increment the counter outside of the conditional
            adaptive_MH_counter += 1;
            burnin_counter += 1;


        }// end of metropolis hastings loop

        // allocate a list to store everything in.
        Rcpp::List ret_list(6);
        ret_list[0] = store_intercepts;
        ret_list[1] = store_coefficients;
        ret_list[2] = store_latent_positions;
        ret_list[3] = store_accept_or_reject;
        ret_list[4] = intercept_proposal_standard_deviations;
        ret_list[5] = LSM_log_likelihood;
        // return everything
        return ret_list;
    }

    // ***********************************************************************//
    //            Sample Data/Parameters From Generative Process              //
    // ***********************************************************************//
    Rcpp::List sample_from_generative_process(arma::vec author_indexes,
                                        arma::cube covariates,
                                        arma::vec alpha_m,
                                        arma::vec beta_n,
                                        bool using_coefficients,
                                        arma::vec intercept_prior_standard_deviations,
                                        arma::vec coefficient_prior_standard_deviations,
                                        arma::vec latent_position_prior_standard_deviations,
                                        int total_number_of_tokens,
                                        int num_documents,
                                        int words_per_doc,
                                        int num_topics,
                                        int num_terms,
                                        int num_actors,
                                        int num_ip,
                                        int num_ld,
                                        arma::vec random_numbers,
                                        Rcpp::List token_topic_assignments,
                                        Rcpp::List token_word_types,
                                        bool resample_word_types,
                                        arma::vec intercepts,
                                        arma::mat coefficients,
                                        arma::cube latent_positions) {

        // sample token topic assignments and word types from generative process
        Rcpp::List ret = sample_token_topics_generative_process(
            token_topic_assignments,
            token_word_types,
            alpha_m,
            beta_n,
            num_documents,
            resample_word_types,
            random_numbers);

        // extract new parameters from the list object
        token_topic_assignments = ret[0];
        token_word_types = ret[1];
        arma::mat document_topic_counts = ret[2];
        arma::vec topic_token_counts = ret[3];
        arma::mat word_type_topic_counts = ret[4];

        // get the interaction pattern parameters. One level up we use the prior
        // means to initialize intercepts, coefficients, latent_positions, at
        // ther mean values, or just pass in their previous values if we are
        // doing backwards sampling. We also need to create vectors of matching
        // standard devaitions in the next level up.
        Rcpp::List ret2 = sample_new_interaction_pattern_parameters(
                intercepts,
                coefficients,
                latent_positions,
                intercept_prior_standard_deviations,
                coefficient_prior_standard_deviations,
                latent_position_prior_standard_deviations,
                using_coefficients);

        // update the input latent positions
        arma::vec temp1 = ret2[0];
        arma::mat temp2 = ret2[1];
        arma::cube temp3 = ret2[2];
        intercepts = temp1;
        coefficients = temp2;
        latent_positions = temp3;

        //generate edge probabilities
        arma::cube edge_probabilities = arma::zeros(num_actors,
                                                    num_actors,
                                                    num_ip);

        // loop over docs, actors, interaction patterns to fill edge_probabilities
        for (int i = 0; i < num_actors; ++i) {
            for (int j = 0; j < num_actors; ++j) {
                // .tube gives us all slices
                arma::vec current_covariates = covariates.tube(i,j);
                if (i != j) {
                    for (int k = 0; k < num_ip; ++k) {
                        edge_probabilities(i,j,k) = mjd::edge_probability(
                            intercepts,
                            coefficients,
                            latent_positions,
                            i,
                            j,
                            current_covariates,
                            k,
                            using_coefficients);
                    }
                }
            }
        }

        // draw interaction patterns from a discrete uniform distribution
        arma::vec interaction_pattern_indices = arma::zeros(num_ip);
        arma::vec interaction_pattern_probs = arma::zeros(num_ip);
        // fill these in with unifrom probabilities
        for (int k = 0; k < num_ip; ++k) {
            interaction_pattern_indices[k] = k;
            interaction_pattern_probs[k] = double(double(1)/double(num_ip));
        }
        //now use RcppArmadillo's implementation of the sample function,
        //with replacement. See:
        //http://gallery.rcpp.org/articles/using-the-Rcpp-based-sample-implementation/
        //for an example and explanation
        arma::vec topic_interaction_patterns = RcppArmadillo::sample(interaction_pattern_indices,
                                                                     num_topics,
                                                                     true,
                                                                     interaction_pattern_probs) ;

        //our last task is to resample edge values
        arma::mat document_edge_matrix = arma::zeros(num_documents,num_actors);
        // loop over documents
        for (int d = 0; d < num_documents; ++d) {
            // we need to allocate row vectors from a matrix before we can
            // convert to arma::vec as desired (typeing issue)
            arma::rowvec temp = document_topic_counts.row(d);
            // now we convert to an arma::vec (ugly but it works)
            arma::vec current_document_topic_counts = arma::conv_to<arma::vec>::from(temp);
            // get the current number of tokens
            int tokens_in_document = arma::sum(current_document_topic_counts);
            int document_sender = author_indexes[d];
            //loop over actors
            for (int r = 0; r < num_actors; ++r) {
                if (r != document_sender) {
                    // now we loop over topics and determine hte edge probability
                    double edge_prob = 0;
                    for (int t = 0; t < num_topics; ++t) {
                        edge_prob += (double(current_document_topic_counts[t])/
                            double(tokens_in_document)) *
                                edge_probabilities(
                                    document_sender,
                                    r,
                                    topic_interaction_patterns[t]);
                    }

                    // now take a random uniform draw and if it is smaller than
                    // the edge prob then we set the edge to 1.
                    double rand_unif_draw = R::runif(0,1);
                    if (edge_prob < rand_unif_draw) {
                        document_edge_matrix(d,r) = 1;
                    }
                }
            }

        }

        // store everything so it can be returned
        Rcpp::List ret_list(10);
        ret_list[0] = token_topic_assignments;
        ret_list[1] = token_word_types;
        ret_list[2] = document_topic_counts;
        ret_list[3] = topic_token_counts;
        ret_list[4] = word_type_topic_counts;
        ret_list[5] = intercepts;
        ret_list[6] = coefficients;
        ret_list[7] = latent_positions;
        ret_list[8] = topic_interaction_patterns;
        ret_list[9] = document_edge_matrix;


        return ret_list;
    }


    // ***********************************************************************//
    //                      Calculate Mean of a Matrix                        //
    // ***********************************************************************//
    double matrix_mean(arma::mat matrix,
                       int subtract_out_entries) {
        int num_cols = matrix.n_cols;
        int num_rows = matrix.n_rows;

        // calculate mean value while optionally holding out some number of
        // entries from the denominator (if they should be zeros, like self ties)
        arma::mat sum_term = arma::sum(matrix);
        double mean_val = double(sum_term(0,0))/double((num_cols*num_rows) - subtract_out_entries);

        return mean_val;
    }

    // ***********************************************************************//
    //                      Calculate Mean of a Vector                        //
    // ***********************************************************************//
    double vector_mean(arma::vec vect) {
        int len = vect.n_elem;

        // calculate mean value
        double mean_val = double(arma::sum(vect))/double(len);

        return mean_val;
    }

    // ***********************************************************************//
    //                      Calculate Sum of LS Distnaces                     //
    // ***********************************************************************//
    double calculate_sum_of_ls_distances(arma::mat positions) {
        int num_actors = positions.n_rows;
        int num_ld = positions.n_cols;
        double sum_term = 0;

        for (int i = 0; i < num_actors; ++i) {
            for (int j = 0; j < num_actors; ++j) {
                if (i < j) {
                    double squared_dist = 0;
                    for (int k = 0; k < num_ld; ++k) {
                        squared_dist += pow((positions(i,k) - positions(j,k)), 2);
                    }
                    sum_term += sqrt(squared_dist);
                }
            }
        }
        return sum_term;
    }

    // ***********************************************************************//
    //         Count Tokens Assigned To Current Interaction Pattern           //
    // ***********************************************************************//
    int count_tokens_assigned_to_ip(arma::vec topic_token_counts,
                                    arma::vec topic_interaction_patterns,
                                    int current_ip) {
        int num_topics = topic_token_counts.n_elem;
        int count = 0;

        for (int i = 0; i < num_topics; ++i) {
            if (topic_interaction_patterns[i] == current_ip) {
                count += topic_token_counts[i];
            }
        }
        return count;
    }

    // ***********************************************************************//
    //            Calculate Statistics for Getting it Right                   //
    // ***********************************************************************//
    arma::vec calculate_statistics_for_getting_it_right(arma::mat document_topic_counts,
        arma::vec topic_token_counts,
        arma::mat word_type_topic_counts,
        arma::vec intercepts,
        arma::mat coefficients,
        arma::cube latent_positions,
        arma::vec topic_interaction_patterns,
        arma::mat document_edge_matrix,
        int num_statistics) {

        // For reference
        //get the number of actors
        //int num_actors = latent_positions.n_cols;
        //get the number of latent dimensions
        //int num_latent_dimensions = latent_positions.n_slices;
        // initialize the proposed values
        //arma::vec proposed_intercepts = arma::zeros(num_interaction_patterns);
        //arma::mat proposed_coefficients = arma::zeros(num_interaction_patterns,
        //                                              num_coefficients);
        //arma::cube proposed_latent_positions = arma::zeros(
        //    num_interaction_patterns,
        //    num_actors,
        //    num_latent_dimensions);

        //variables we will need for calcualting statistics
        int number_of_interaction_patterns = intercepts.n_elem;
        int number_of_topics = topic_token_counts.n_elem;
        int number_of_word_types = word_type_topic_counts.n_rows;
        int num_docs = document_topic_counts.n_rows;
        int num_actors = latent_positions.n_cols;

        // allocate a vector in which to store statistics we calculate on our
        // data and latent variables
        arma::vec statistics = arma::zeros(num_statistics);
        int stat_counter = 0;

        // For reference, this is what we are saving.
        //we want average ls positions, coefficients, the intercept + sum of
        //distances for each interaction pattern, number of tokens assigned to each
        //interaction pattern
        //number_of_statistics += 5*num_ip;
        //number_of_statistics += num_terms + num_topics;
        //number_of_statistics += 2; //average cluster assignment and mean network

        // all interaction pattern specific stats

        // intercept for each interaction pattern
        for (int i = 0; i < number_of_interaction_patterns; ++i) {
            statistics[stat_counter] = intercepts[i];
            stat_counter += 1;
        }

        // mean of coefficients positions in each interaction pattern
        for (int i = 0; i < number_of_interaction_patterns; ++i) {
            arma::vec temp = coefficients.row(i);
            statistics[stat_counter] = vector_mean(temp);
            stat_counter += 1;
        }

        // mean of LS positions in each interaction pattern
        for (int i = 0; i < number_of_interaction_patterns; ++i) {
            arma::mat temp = latent_positions.tube(i,0,i,num_actors-1);
            statistics[stat_counter] = matrix_mean(temp,
                                                   0);
            stat_counter += 1;
        }

        // sum of LS distances in each interaction pattern
        for (int i = 0; i < number_of_interaction_patterns; ++i) {
            arma::mat temp = latent_positions.tube(i,0,i,num_actors-1);
            statistics[stat_counter] = calculate_sum_of_ls_distances(temp);
            stat_counter += 1;
        }

        // sum of tokens assigned to topics in each interaction pattern
        for (int i = 0; i < number_of_interaction_patterns; ++i) {
            statistics[stat_counter] = count_tokens_assigned_to_ip(topic_token_counts,
                                                                   topic_interaction_patterns,
                                                                   i);
            stat_counter += 1;
        }

        // all topic specific stats
        // number of tokens assigned to each topic
        for (int i = 0; i < number_of_topics; ++i) {
            statistics[stat_counter] = topic_token_counts[i];
            stat_counter += 1;
        }

        // all word-type specific stats
        // number of tokens assigned to each word type
        for (int i = 0; i < number_of_word_types; ++i) {
            arma::vec temp = word_type_topic_counts.row(i);
            statistics[stat_counter] = arma::sum(temp);
            stat_counter += 1;
        }

        // now get density which will be calcualted as the fraction of all
        // possible email recipients
        statistics[stat_counter] = matrix_mean(document_edge_matrix,
                                               num_docs);
        stat_counter += 1;

        // now get the mean of cluster assignments
        statistics[stat_counter] = vector_mean(topic_interaction_patterns);
        stat_counter += 1;

        return statistics;
    }

} // end of MJD namespace


    // ***********************************************************************//
    // ***********************************************************************//
    // ***********************************************************************//
    //                                TEST FUNCTIONS:                         //
    //  These functions expose the c++ functions to R so that we can test them//
    //  easily.                                                               //
    // ***********************************************************************//
    // ***********************************************************************//
    // ***********************************************************************//







using namespace Rcpp;
// for testing we will wrap and export this function so it is available in R.
// [[Rcpp::export]]
double ep(arma::vec intercepts,
       arma::mat coefficients,
       arma::cube latent_positions,
       int sender,
       int recipient,
       arma::vec current_covariates,
       int interaction_pattern,
       bool using_coefficients){

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
int lsms(arma::vec unnormalized_discrete_distribution, int seed, double u=NA_REAL){

    // Set RNG and define uniform distribution
    boost::mt19937 generator(seed);
    boost::uniform_01<double> uniform_distribution;

    // get the random uniform draw
    double rand_num = uniform_distribution(generator);

    if (arma::is_finite(u)) {
      rand_num = u;
    }

    // take a draw from the unnormalized log distribution
    int temp = mjd::log_space_multinomial_sampler (
            unnormalized_discrete_distribution,
            rand_num);
    //increment by one to go back to 1 indexing in R.
    return (temp + 1);
}

// [[Rcpp::export]]
double sotep(arma::cube edge_probabilities,
          int tokens_in_document,
          int current_token_topic_assignment,
          arma::vec current_document_topic_counts,
          bool leave_out_current_token,
          arma::vec topic_interaction_patterns,
          int document_sender,
          int document_recipient,
          int leave_out_topic){

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
double ppipp(arma::vec intercepts,
          arma::mat coefficients,
          arma::cube latent_positions,
          double intercept_prior_mean,
          double intercept_prior_standard_deviation,
          double coefficient_prior_mean,
          double coefficient_prior_standard_deviation,
          double latent_position_prior_mean,
          double latent_position_prior_standard_deviation,
          bool using_coefficients){

    double log_prob = mjd::prior_probability_interaction_pattern_parameters (
        intercepts,
        coefficients,
        latent_positions,
        intercept_prior_mean,
        intercept_prior_standard_deviation,
        coefficient_prior_mean,
        coefficient_prior_standard_deviation,
        latent_position_prior_mean,
        latent_position_prior_standard_deviation,
        using_coefficients);

    return log_prob;
}

// [[Rcpp::export]]
List snipp(arma::vec intercepts,
          arma::mat coefficients,
          arma::cube latent_positions,
          arma::vec intercept_proposal_standard_deviations,
          arma::vec coefficient_proposal_standard_deviations,
          arma::vec latent_position_proposal_standard_deviations,
          bool using_coefficients){

    List new_params = mjd::sample_new_interaction_pattern_parameters (
        intercepts,
        coefficients,
        latent_positions,
        intercept_proposal_standard_deviations,
        coefficient_proposal_standard_deviations,
        latent_position_proposal_standard_deviations,
        using_coefficients);

    return new_params;
}

// [[Rcpp::export]]
double lsmc(arma::cube edge_probabilities,
            int tokens_in_document,
            int topic,
            int current_token_topic_assignment,
            arma::vec current_document_topic_counts,
            arma::vec document_edge_values,
            arma::vec topic_interaction_patterns,
            int document_sender){

    double contrib = mjd::lsm_contribution (
        edge_probabilities,
        tokens_in_document,
        topic,
        current_token_topic_assignment,
        current_document_topic_counts,
        document_edge_values,
        topic_interaction_patterns,
        document_sender);

    return contrib;
}



// [[Rcpp::export]]
double ldac(int tokens_in_document,
          int current_token_topic_assignment,
          arma::vec current_document_topic_counts,
          arma::mat word_type_topic_counts,
          arma::vec topic_token_counts,
          int topic,
          int current_word_type,
          arma::vec alpha_m,
          arma::vec beta_n,
          double beta){

    double contribution = mjd::lda_contribution(
        tokens_in_document,
        current_token_topic_assignment,
        current_document_topic_counts,
        word_type_topic_counts,
        topic_token_counts,
        topic,
        current_word_type,
        alpha_m,
        beta_n,
        beta);

    return contribution;
}


// [[Rcpp::export]]
arma::vec ustta(arma::cube edge_probabilities,
        int tokens_in_document,
        int current_token_topic_assignment,
        arma::vec current_document_topic_counts,
        arma::mat word_type_topic_counts,
        arma::vec topic_token_counts,
        int current_word_type,
        arma::vec alpha_m,
        arma::vec beta_n,
        arma::vec document_edge_values,
        arma::vec topic_interaction_patterns,
        int document_sender,
        double rand_num,
        bool parallel,
        bool use_cached_token_topic_distribution,
        arma::vec cached_token_topic_distribution){

   arma::vec assignment = mjd::update_single_token_topic_assignment(
            edge_probabilities,
            tokens_in_document,
            current_token_topic_assignment,
            current_document_topic_counts,
            word_type_topic_counts,
            topic_token_counts,
            current_word_type,
            alpha_m,
            beta_n,
            document_edge_values,
            topic_interaction_patterns,
            document_sender,
            rand_num,
            parallel,
            use_cached_token_topic_distribution,
            cached_token_topic_distribution);

    return assignment;
}

// [[Rcpp::export]]
List utta(arma::vec author_indexes,
        arma::mat document_edge_matrix,
        arma::vec topic_interaction_patterns,
        arma::mat document_topic_counts,
        arma::mat word_type_topic_counts,
        arma::vec topic_token_counts,
        Rcpp::List token_topic_assignments,
        Rcpp::List token_word_types,
        arma::vec intercepts,
        arma::mat coefficients,
        arma::cube latent_positions,
        arma::cube covariates,
        arma::vec alpha_m,
        arma::vec beta_n,
        arma::vec random_numbers,
        bool using_coefficients,
        bool parallel){

    //random_numbers has length equal to the total number of tokens in the corpus.
    List ret_list = mjd::update_token_topic_assignments(
            author_indexes,
            document_edge_matrix,
            topic_interaction_patterns,
            document_topic_counts,
            word_type_topic_counts,
            topic_token_counts,
            token_topic_assignments,
            token_word_types,
            intercepts,
            coefficients,
            latent_positions,
            covariates,
            alpha_m,
            beta_n,
            random_numbers,
            using_coefficients,
            parallel);

    return ret_list;

}


// [[Rcpp::export]]
List uipp(arma::vec author_indexes,
          arma::mat document_edge_matrix,
          arma::mat document_topic_counts,
          arma::vec topic_interaction_patterns,
          arma::vec intercepts,
          arma::mat coefficients,
          arma::cube latent_positions,
          arma::cube covariates,
          bool using_coefficients,
          double intercept_prior_mean,
          double intercept_prior_standard_deviation,
          arma::vec intercept_proposal_standard_deviations,
          double coefficient_prior_mean,
          double coefficient_prior_standard_deviation,
          arma::vec coefficient_proposal_standard_deviations,
          double latent_position_prior_mean,
          double latent_position_prior_standard_deviation,
          arma::vec latent_position_proposal_standard_deviations,
          double random_number,
          arma::cube edge_probabilities){

    //random_numbers has length equal to the total number of tokens in the corpus.
    List ret_list =  mjd::update_interaction_pattern_parameters(
            author_indexes,
            document_edge_matrix,
            document_topic_counts,
            topic_interaction_patterns,
            intercepts,
            coefficients,
            latent_positions,
            covariates,
            using_coefficients,
            intercept_prior_mean,
            intercept_prior_standard_deviation,
            intercept_proposal_standard_deviations,
            coefficient_prior_mean,
            coefficient_prior_standard_deviation,
            coefficient_proposal_standard_deviations,
            latent_position_prior_mean,
            latent_position_prior_standard_deviation,
            latent_position_proposal_standard_deviations,
            random_number,
            edge_probabilities);

    return ret_list;

}


// [[Rcpp::export]]
arma::vec utipa(arma::vec author_indexes,
           arma::mat document_edge_matrix,
           arma::mat document_topic_counts,
           arma::vec topic_interaction_patterns,
           arma::vec intercepts,
           arma::mat coefficients,
           arma::cube latent_positions,
           arma::cube covariates,
           bool using_coefficients,
           arma::vec random_numbers,
           arma::cube edge_probabilities){

    //random_numbers has length equal to the total number of tokens in the corpus.
    arma::vec return_vec =  mjd::update_topic_interaction_pattern_assignments(
            author_indexes,
            document_edge_matrix,
            document_topic_counts,
            topic_interaction_patterns,
            intercepts,
            coefficients,
            latent_positions,
            covariates,
            using_coefficients,
            random_numbers,
            edge_probabilities);

    return return_vec;

}


// [[Rcpp::export]]
List am(arma::vec intercept_proposal_standard_deviations,
        arma::vec coefficient_proposal_standard_deviations,
        arma::vec latent_position_proposal_standard_deviations,
        arma::vec accept_rates,
        double target_accept_rate,
        double tollerance,
        double update_size){

    List ret_list = mjd::adaptive_metropolis(
            intercept_proposal_standard_deviations,
            coefficient_proposal_standard_deviations,
            latent_position_proposal_standard_deviations,
            accept_rates,
            target_accept_rate,
            tollerance,
            update_size);

    return ret_list;

}

// [[Rcpp::export]]
List model_inference(arma::vec author_indexes,
                     arma::mat document_edge_matrix,
                     arma::mat document_topic_counts,
                     arma::vec topic_interaction_patterns,
                     arma::mat word_type_topic_counts,
                     arma::vec topic_token_counts,
                     Rcpp::List token_topic_assignments,
                     Rcpp::List token_word_types,
                     arma::vec intercepts,
                     arma::mat coefficients,
                     arma::cube latent_positions,
                     arma::cube covariates,
                     arma::vec alpha_m,
                     arma::vec beta_n,
                     bool using_coefficients,
                     double intercept_prior_mean,
                     double intercept_prior_standard_deviation,
                     arma::vec intercept_proposal_standard_deviations,
                     double coefficient_prior_mean,
                     double coefficient_prior_standard_deviation,
                     arma::vec coefficient_proposal_standard_deviations,
                     double latent_position_prior_mean,
                     double latent_position_prior_standard_deviation,
                     arma::vec latent_position_proposal_standard_deviations,
                     double target_accept_rate,
                     double tollerance,
                     double update_size,
                     int seed,
                     int iterations,
                     int metropolis_iterations,
                     int total_number_of_tokens,
                     int iterations_before_t_i_p_updates,
                     int update_t_i_p_every_x_iterations,
                     bool perform_adaptive_metropolis,
                     int slice_sample_every_x_iterations,
                     double slice_sample_step_size,
                     bool parallel){

    List ret_list =  mjd::inference(
        author_indexes,
        document_edge_matrix,
        document_topic_counts,
        topic_interaction_patterns,
        word_type_topic_counts,
        topic_token_counts,
        token_topic_assignments,
        token_word_types,
        intercepts,
        coefficients,
        latent_positions,
        covariates,
        alpha_m,
        beta_n,
        using_coefficients,
        intercept_prior_mean,
        intercept_prior_standard_deviation,
        intercept_proposal_standard_deviations,
        coefficient_prior_mean,
        coefficient_prior_standard_deviation,
        coefficient_proposal_standard_deviations,
        latent_position_prior_mean,
        latent_position_prior_standard_deviation,
        latent_position_proposal_standard_deviations,
        target_accept_rate,
        tollerance,
        update_size,
        seed,
        iterations,
        metropolis_iterations,
        total_number_of_tokens,
        iterations_before_t_i_p_updates,
        update_t_i_p_every_x_iterations,
        perform_adaptive_metropolis,
        slice_sample_every_x_iterations,
        slice_sample_step_size,
        parallel);

    return ret_list;

}

// [[Rcpp::export]]
arma::vec mjd_rdirichlet(arma::vec alpha_m) {
    arma::vec to_return = mjd::rdirichlet(alpha_m);
    return(to_return);
}

// [[Rcpp::export]]
Rcpp::List sttgp(
        Rcpp::List token_topic_assignments,
        Rcpp::List token_word_types,
        arma::vec alpha_m,
        arma::vec beta_n,
        int number_of_documents,
        bool resample_word_types,
        arma::vec random_numbers) {

    List ret_list = mjd::sample_token_topics_generative_process(
            token_topic_assignments,
            token_word_types,
            alpha_m,
            beta_n,
            number_of_documents,
            resample_word_types,
            random_numbers);

    return ret_list;
}

// [[Rcpp::export]]
Rcpp::List mh_to_convergence(
        arma::vec author_indexes,
        arma::mat document_edge_matrix,
        arma::mat document_topic_counts,
        arma::vec topic_interaction_patterns,
        arma::vec intercepts,
        arma::mat coefficients,
        arma::cube latent_positions,
        arma::cube covariates,
        bool using_coefficients,
        double intercept_prior_mean,
        double intercept_prior_standard_deviation,
        arma::vec intercept_proposal_standard_deviations,
        double coefficient_prior_mean,
        double coefficient_prior_standard_deviation,
        arma::vec coefficient_proposal_standard_deviations,
        double latent_position_prior_mean,
        double latent_position_prior_standard_deviation,
        arma::vec latent_position_proposal_standard_deviations,
        double target_accept_rate,
        double tollerance,
        double update_size,
        int seed,
        int metropolis_iterations,
        int adaptive_metropolis_every_x_iterations,
        int stop_adaptive_metropolis_after_x_updates,
        int samples_to_store,
        int sample_every,
        int burnin) {

    List ret_list = mjd::run_metropolis_hastings_to_convergence(
        author_indexes,
        document_edge_matrix,
        document_topic_counts,
        topic_interaction_patterns,
        intercepts,
        coefficients,
        latent_positions,
        covariates,
        using_coefficients,
        intercept_prior_mean,
        intercept_prior_standard_deviation,
        intercept_proposal_standard_deviations,
        coefficient_prior_mean,
        coefficient_prior_standard_deviation,
        coefficient_proposal_standard_deviations,
        latent_position_prior_mean,
        latent_position_prior_standard_deviation,
        latent_position_proposal_standard_deviations,
        target_accept_rate,
        tollerance,
        update_size,
        seed,
        metropolis_iterations,
        adaptive_metropolis_every_x_iterations,
        stop_adaptive_metropolis_after_x_updates,
        samples_to_store,
        sample_every,
        burnin);

    return ret_list;
}


// [[Rcpp::export]]
arma::mat gir(arma::vec author_indexes,
         arma::cube covariates,
         arma::vec alpha_m,
         arma::vec beta_n,
         bool using_coefficients,
         double intercept_prior_mean,
         double intercept_prior_standard_deviation,
         arma::vec intercept_proposal_standard_deviations,
         double coefficient_prior_mean,
         double coefficient_prior_standard_deviation,
         arma::vec coefficient_proposal_standard_deviations,
         double latent_position_prior_mean,
         double latent_position_prior_standard_deviation,
         arma::vec latent_position_proposal_standard_deviations,
         double target_accept_rate,
         double tollerance,
         double update_size,
         int seed,
         int iterations,
         int metropolis_iterations,
         int total_number_of_tokens,
         int iterations_before_t_i_p_updates,
         int update_t_i_p_every_x_iterations,
         bool perform_adaptive_metropolis,
         int slice_sample_every_x_iterations,
         double slice_sample_step_size,
         bool parallel,
         int num_documents,
         int words_per_doc,
         int num_topics,
         int num_terms,
         int num_actors,
         int num_ip,
         int num_ld,
         int GiR_samples,
         bool forward_sample,
         Rcpp::List token_topic_assignments,
         Rcpp::List token_word_types,
         bool resample_word_types){

    // Set RNG and define uniform distribution
    boost::mt19937 generator(seed);
    boost::uniform_01<double> uniform_distribution;

    // we will calculate all statistics on which we wish to compare the two
    // chains in C++ and will then store them in a matrix.
    // We pre-allocate it now, and then return them in a list object.
    int number_of_statistics = 0;

    //we want average ls positions, coefficients, the intercept + sum of
    //distances for each interaction pattern, number of tokens assigned to each
    //interaction pattern
    number_of_statistics += 5*num_ip;
    number_of_statistics += num_terms + num_topics;
    number_of_statistics += 2; //average cluster assignment and mean network density.

    arma::mat sample_statistics = arma::zeros(GiR_samples,number_of_statistics);

    // determine hte number of covariates we are using so we can determine
    // the correct number of coefficients
    int num_coefficients = 2;
    if(using_coefficients){
        num_coefficients = covariates.n_slices;
    }
    // allocate interaction pattern parameter and prior sd vectors.
    arma::vec intercept_prior_standard_deviations = arma::zeros(num_ip);
    arma::vec coefficient_prior_standard_deviations = arma::zeros(num_ip);
    arma::vec latent_position_prior_standard_deviations = arma::zeros(num_ip);
    arma::vec intercepts = arma::zeros(num_ip);
    arma::mat coefficients = arma::zeros(num_ip,num_coefficients);
    arma::cube latent_positions = arma::zeros(num_actors,num_actors,num_ld);
    //now fill them in with the appropriate values
    intercept_prior_standard_deviations.fill(intercept_prior_standard_deviation);
    coefficient_prior_standard_deviations.fill(coefficient_prior_standard_deviation);
    latent_position_prior_standard_deviations.fill(latent_position_prior_standard_deviation);
    intercepts.fill(intercept_prior_mean);
    coefficients.fill(coefficient_prior_mean);
    latent_positions.fill(latent_position_prior_mean);

    // allocate data structures to store
    if (forward_sample) {
        // Forward Samples:
        for (int i = 0; i < GiR_samples; ++i) {
            // Take a draw from the generative process
            // generate a vector of random numbers to pass in to the topic-token
            // update function.
            arma::vec random_numbers = arma::zeros(total_number_of_tokens);
            for (int k = 0; k < total_number_of_tokens; ++k) {
                random_numbers[k] = uniform_distribution(generator);
            }

            Rcpp::List ret = mjd::sample_from_generative_process(author_indexes,
                  covariates,
                  alpha_m,
                  beta_n,
                  using_coefficients,
                  intercept_prior_standard_deviations,
                  coefficient_prior_standard_deviations,
                  latent_position_prior_standard_deviations,
                  total_number_of_tokens,
                  num_documents,
                  words_per_doc,
                  num_topics,
                  num_terms,
                  num_actors,
                  num_ip,
                  num_ld,
                  random_numbers,
                  token_topic_assignments,
                  token_word_types,
                  resample_word_types,
                  intercepts,
                  coefficients,
                  latent_positions);

            //now extract everything from the list.
            Rcpp::List token_topic_assignments = ret[0];
            Rcpp::List token_word_types = ret[1];
            arma::mat document_topic_counts = ret[2];
            arma::vec topic_token_counts = ret[3];
            arma::mat word_type_topic_counts = ret[4];
            arma::vec temp1 = ret[5];
            intercepts = temp1;
            arma::mat temp2 = ret[6];
            coefficients = temp2;
            arma::cube temp3 = ret[7];
            latent_positions= temp3;
            arma::vec topic_interaction_patterns = ret[8];
            arma::mat document_edge_matrix = ret[9];

            // Calculate statistics
            arma::vec stats =  mjd::calculate_statistics_for_getting_it_right(document_topic_counts,
                topic_token_counts,
                word_type_topic_counts,
                intercepts,
                coefficients,
                latent_positions,
                topic_interaction_patterns,
                document_edge_matrix,
                number_of_statistics);

            // Store statistics
            for (int k = 0; k < number_of_statistics; ++k) {
                sample_statistics(i,k) = stats[k];
            }

        }
    } else {
        // Backward Samples:
        // Run the generative process for the first time and extract relevant
        // datastructures
        // Take a draw from the generative process
        // generate a vector of random numbers to pass in to the topic-token
        // update function.
        arma::vec random_numbers = arma::zeros(total_number_of_tokens);
        for (int k = 0; k < total_number_of_tokens; ++k) {
            random_numbers[k] = uniform_distribution(generator);
        }

        Rcpp::List ret = mjd::sample_from_generative_process(author_indexes,
                                                             covariates,
                                                             alpha_m,
                                                             beta_n,
                                                             using_coefficients,
                                                             intercept_prior_standard_deviations,
                                                             coefficient_prior_standard_deviations,
                                                             latent_position_prior_standard_deviations,
                                                             total_number_of_tokens,
                                                             num_documents,
                                                             words_per_doc,
                                                             num_topics,
                                                             num_terms,
                                                             num_actors,
                                                             num_ip,
                                                             num_ld,
                                                             random_numbers,
                                                             token_topic_assignments,
                                                             token_word_types,
                                                             resample_word_types,
                                                             intercepts,
                                                             coefficients,
                                                             latent_positions);

        //now extract everything from the list.
        Rcpp::List token_topic_assignments = ret[0];
        Rcpp::List token_word_types = ret[1];
        arma::mat document_topic_counts = ret[2];
        arma::vec topic_token_counts = ret[3];
        arma::mat word_type_topic_counts = ret[4];
        arma::vec temp1 = ret[5];
        intercepts = temp1;
        arma::mat temp2 = ret[6];
        coefficients = temp2;
        arma::cube temp3 = ret[7];
        latent_positions= temp3;
        arma::vec topic_interaction_patterns = ret[8];
        arma::mat document_edge_matrix = ret[9];

        for (int i = 0; i < GiR_samples; ++i) {
            // Run inference for 5 Gibbs (50 MH per Gibbs) iterations.

            Rcpp::List ret_list =  mjd::inference(
                author_indexes,
                document_edge_matrix,
                document_topic_counts,
                topic_interaction_patterns,
                word_type_topic_counts,
                topic_token_counts,
                token_topic_assignments,
                token_word_types,
                intercepts,
                coefficients,
                latent_positions,
                covariates,
                alpha_m,
                beta_n,
                using_coefficients,
                intercept_prior_mean,
                intercept_prior_standard_deviation,
                intercept_proposal_standard_deviations,
                coefficient_prior_mean,
                coefficient_prior_standard_deviation,
                coefficient_proposal_standard_deviations,
                latent_position_prior_mean,
                latent_position_prior_standard_deviation,
                latent_position_proposal_standard_deviations,
                target_accept_rate,
                tollerance,
                update_size,
                seed,
                iterations,
                metropolis_iterations,
                total_number_of_tokens,
                iterations_before_t_i_p_updates,
                update_t_i_p_every_x_iterations,
                perform_adaptive_metropolis,
                slice_sample_every_x_iterations,
                slice_sample_step_size,
                parallel);

            //extract parameters
            arma::mat inf1 = ret_list[4];
            document_topic_counts = inf1;
            arma::mat inf2 = ret_list[5];
            word_type_topic_counts = inf2;
            arma::mat inf3 = ret_list[6];
            topic_token_counts = inf3;
            arma::mat inf4 = ret_list[7];
            token_topic_assignments = inf4;
            arma::vec inf5 = ret_list[11];
            intercepts = inf5;
            arma::mat inf6 = ret_list[12];
            coefficients = inf6;
            arma::cube inf7 = ret_list[13];
            latent_positions= inf7;
            arma::vec inf8 = ret_list[14];
            topic_interaction_patterns = inf8;
            arma::mat inf9 = ret_list[15];
            word_type_topic_counts = inf9;

            // Calculate statistics
            arma::vec stats =  mjd::calculate_statistics_for_getting_it_right(document_topic_counts,
                                                                              topic_token_counts,
                                                                              word_type_topic_counts,
                                                                              intercepts,
                                                                              coefficients,
                                                                              latent_positions,
                                                                              topic_interaction_patterns,
                                                                              document_edge_matrix,
                                                                              number_of_statistics);

            // Store statistics
            for (int k = 0; k < number_of_statistics; ++k) {
                sample_statistics(i,k) = stats[k];
            }

            // Take a draw from the generative process using updated parameters.
            // generate a vector of random numbers to pass in to the topic-token
            // update function.
            arma::vec random_numbers = arma::zeros(total_number_of_tokens);
            for (int k = 0; k < total_number_of_tokens; ++k) {
                random_numbers[k] = uniform_distribution(generator);
            }

            Rcpp::List ret = mjd::sample_from_generative_process(author_indexes,
                                                                 covariates,
                                                                 alpha_m,
                                                                 beta_n,
                                                                 using_coefficients,
                                                                 intercept_prior_standard_deviations,
                                                                 coefficient_prior_standard_deviations,
                                                                 latent_position_prior_standard_deviations,
                                                                 total_number_of_tokens,
                                                                 num_documents,
                                                                 words_per_doc,
                                                                 num_topics,
                                                                 num_terms,
                                                                 num_actors,
                                                                 num_ip,
                                                                 num_ld,
                                                                 random_numbers,
                                                                 token_topic_assignments,
                                                                 token_word_types,
                                                                 resample_word_types,
                                                                 intercepts,
                                                                 coefficients,
                                                                 latent_positions);

            //now extract everything from the list.
            Rcpp::List temp4 = ret[0];
            token_topic_assignments = temp4;
            Rcpp::List temp5 = ret[1];
            token_word_types = temp5;
            arma::mat temp6 = ret[2];
            document_topic_counts = temp6;
            arma::vec temp7 = ret[3];
            topic_token_counts = temp7;
            arma::mat temp8 = ret[4];
            word_type_topic_counts = temp8;
            arma::vec temp11 = ret[5];
            intercepts = temp11;
            arma::mat temp12 = ret[6];
            coefficients = temp12;
            arma::cube temp13 = ret[7];
            latent_positions= temp13;
            arma::vec temp14 = ret[8];
            topic_interaction_patterns = temp14;
            arma::mat temp15 = ret[9];
            document_edge_matrix = temp15;

        }


    }

    // return the statistics
    return sample_statistics;
}


