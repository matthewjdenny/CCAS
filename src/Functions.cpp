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
        for (int i = 0; i < num_latent_spaces; ++i) {
            linear_predictor += abs((latent_positions(interaction_pattern, sender, i) -
                latent_positions(interaction_pattern, recipient, i)));
        }

        // now run through logit
        double edge_prob = (1/double(1+ exp(-linear_predictor)));
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

	// ***********************************************************************//
	//         Get Prior Probability of Interaction Pattern Parameters        //
	// ***********************************************************************//

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

    // ***********************************************************************//
    //                Sample New Interaction Pattern Parameters               //
    // ***********************************************************************//

    Rcpp::List sample_new_interaction_pattern_parameters(
        arma::vec intercepts,
        arma::mat coefficients,
        arma::cube latent_positions,
        arma::vec intercept_proposal_variances,
        arma::vec coefficient_proposal_variances,
        arma::vec latent_position_proposal_variances,
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
                                              intercept_proposal_variances[i]);
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
                                                    coefficient_proposal_variances[i]);
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
                        latent_position_proposal_variances[i]);
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

        // get the total count of tokens currently assocaited with topic
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
    //                 Update Single Token Topic Assignment                   //
    // ***********************************************************************//

    int update_single_token_topic_assignment(
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
            double rand_num) {

        // calculate beta, the sum over beta_n
        double beta = arma::sum(beta_n);

        // get the number of topics
        int number_of_topics = alpha_m.n_elem;

        // generate a zero vector of length (number of topics)
        arma::vec unnormalized_log_distribution = arma::zeros(number_of_topics);

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

            // add the values (since we are working in log space) and put them
            // in the appropriate bin the in the distribution.
            unnormalized_log_distribution[i] = lsm_contr + lda_contr;
        }

        // now we need to sample from this distribution
        int new_assignment = mjd::log_space_multinomial_sampler (
            unnormalized_log_distribution,
            rand_num);

        return new_assignment;
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
            bool using_coefficients) {

        // get important constants
        int number_of_documents = document_edge_matrix.n_rows;
        int number_of_actors = document_edge_matrix.n_cols;
        int number_of_interaction_patterns = intercepts.n_elem;
        int rand_num_counter = 0;

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
            // get the current token topic assignments as a vector
            arma::vec current_token_topic_assignments = token_topic_assignments[i];
            // get the current token word types as a vector
            arma::vec current_token_word_types = token_word_types[i];
            // get the current number of tokens
            int tokens_in_document = current_token_topic_assignments.n_elem;
            // allocate all of our document specific variables
            arma::vec current_document_topic_counts = document_topic_counts.row(i);
            arma::vec document_edge_values = document_edge_matrix.row(i);
            int document_sender = author_indexes[i];

            // loop over tokens
            for (int j = 0; j < tokens_in_document; ++j) {
                // allocate all of our token specific variables
                int current_token_topic_assignment = current_token_topic_assignments[j];
                int current_word_type = current_token_word_types[j];
                double rand_num = random_numbers[rand_num_counter];
                rand_num_counter += 1;

                // now get the new assignment
                int new_topic_assignment = update_single_token_topic_assignment(
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
                    rand_num);

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
            double intercept_prior_variance,
            arma::vec intercept_proposal_variances,
            double coefficient_prior_mean,
            double coefficient_prior_variance,
            arma::vec coefficient_proposal_variances,
            double latent_position_prior_mean,
            double latent_position_prior_variance,
            arma::vec latent_position_proposal_variances,
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
            intercept_proposal_variances,
            coefficient_proposal_variances,
            latent_position_proposal_variances,
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
            intercept_prior_variance,
            coefficient_prior_mean,
            coefficient_prior_variance,
            latent_position_prior_mean,
            latent_position_prior_variance,
            using_coefficients);

        // get the log prior probability of the proposed interaction pattern
        // parameters.
        double log_proposed_prior = mjd::prior_probability_interaction_pattern_parameters(
            proposed_intercepts,
            proposed_coefficients,
            proposed_latent_positions,
            intercept_prior_mean,
            intercept_prior_variance,
            coefficient_prior_mean,
            coefficient_prior_variance,
            latent_position_prior_mean,
            latent_position_prior_variance,
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
            arma::vec current_document_topic_counts = document_topic_counts.row(i);
            arma::vec document_edge_values = document_edge_matrix.row(i);
            // get the current number of tokens
            int tokens_in_document = arma::sum(current_document_topic_counts);
            int document_sender = author_indexes[i];
            // loop over tokens
            for (int j = 0; j < number_of_actors; ++j) {
                // if the assignment changed, then we need to update everything.
                if (document_sender != j) {
                    if (document_edge_values[i] == 1) {

                        double temp = sum_over_t_edge_probability (
                            edge_probabilities,
                            tokens_in_document,
                            -1,
                            current_document_topic_counts,
                            false,
                            topic_interaction_patterns,
                            document_sender,
                            j,
                            -1);
                        log_current_probability += log(temp);

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
                        log_proposed_probability += log(temp2);
                    } else {
                        double temp = sum_over_t_edge_probability (
                            edge_probabilities,
                            tokens_in_document,
                            -1,
                            current_document_topic_counts,
                            false,
                            topic_interaction_patterns,
                            document_sender,
                            j,
                            -1);
                        log_current_probability += log(1- temp);
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
                        log_proposed_probability += log(1- temp2);
                    }
                } // end of condition making sure actor is not author
            } // end of loop over actors
        }// end of loop over documents

        double accept_log_prob = log_proposed_probability + log_proposed_prior -
             log_current_probability - log_current_prior;
        double log_random_number = log(random_number);
        Rcpp::List to_return(5);

        if (accept_log_prob > log_random_number) {
            to_return[0] = proposed_intercepts;
            to_return[1] = proposed_coefficients;
            to_return[2] = proposed_latent_positions;
            to_return[3] = proposed_edge_probabilities;
            to_return[4] = 1; //tells us whether we accepted
        } else {
            to_return[0] = intercepts;
            to_return[1] = coefficients;
            to_return[2] = latent_positions;
            to_return[3] = edge_probabilities;
            to_return[4] = 0; //tells us whether we accepted
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
                arma::vec current_document_topic_counts = document_topic_counts.row(i);
                arma::vec document_edge_values = document_edge_matrix.row(i);
                // get the current number of tokens
                int tokens_in_document = arma::sum(current_document_topic_counts);
                int document_sender = author_indexes[i];
                // loop over tokens
                for (int j = 0; j < number_of_actors; ++j) {
                    // if the assignment changed, then we need to update everything.
                    if (document_sender != j) {
                        if (document_edge_values[i] == 1) {
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
                            held_out_sum_over_t_terms(i,j) = temp;
                        } else {
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
                            held_out_sum_over_t_terms(i,j) = 1- temp;
                        }
                    } // end of condition making sure actor is not author
                } // end of loop over actors
            }// end of loop over documents

            // loop over clusters to populate distribution
            for (int c = 0; c < number_of_interaction_patterns; ++c) {
                double log_prob = 0;
                // loop over documents
                for (int i = 0; i < number_of_documents; ++i) {
                    // allocate all of our document specific variables
                    arma::vec current_document_topic_counts = document_topic_counts.row(i);
                    arma::vec document_edge_values = document_edge_matrix.row(i);
                    // get the current number of tokens
                    int tokens_in_document = arma::sum(current_document_topic_counts);
                    int document_sender = author_indexes[i];
                    // loop over tokens
                    for (int j = 0; j < number_of_actors; ++j) {
                        // if the assignment changed, then we need to update everything.
                        if (document_sender != j) {
                            double ttc = double(current_document_topic_counts[t])/
                            double(tokens_in_document);
                            if (document_edge_values[i] == 1) {
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

    // get the random uniform draw
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
          arma::vec intercept_proposal_variances,
          arma::vec coefficient_proposal_variances,
          arma::vec latent_position_proposal_variances,
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
        intercept_proposal_variances,
        coefficient_proposal_variances,
        latent_position_proposal_variances,
        using_coefficients);

    return new_params;
}

// [[Rcpp::export]]
double lsmc(NumericVector edge_probs,
            int tokens_in_document,
            int topic,
            int current_token_topic_assignment,
            arma::vec current_document_topic_counts,
            arma::vec document_edge_values,
            arma::vec topic_interaction_patterns,
            int document_sender){

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
int ustta(NumericVector edge_probs,
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
        double rand_num){

    // we have to do this stupid trick to pass in 3d arrays from R. We pass in as
    // a vector, then instatiate a cube object from there.
    IntegerVector arrayDims = edge_probs.attr("dim");
    arma::cube edge_probabilities(edge_probs.begin(), arrayDims[0], arrayDims[1],
                                  arrayDims[2], false);

    int assignment = mjd::update_single_token_topic_assignment(
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
            rand_num);

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
        NumericVector latent_pos,
        NumericVector covars,
        arma::vec alpha_m,
        arma::vec beta_n,
        arma::vec random_numbers,
        bool using_coefficients){

    // Make sure you supply a NumericVector as input!
    // we have to do this stupid trick to pass in 3d arrays from R. We pass in as
    // a vector, then instatiate a cube object from there.
    IntegerVector arrayDims = latent_pos.attr("dim");
    arma::cube latent_positions(latent_pos.begin(), arrayDims[0], arrayDims[1], arrayDims[2], false);

    IntegerVector arrayDims2 = covars.attr("dim");
    arma::cube covariates(covars.begin(), arrayDims2[0], arrayDims2[1], arrayDims2[2], false);

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
            using_coefficients);

    return ret_list;

}





