// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(BH)]]

#include <RcppArmadillo.h>

using std::log;
using std::exp;
using std::max;
using std::abs;

namespace mjd {
    // this is super annoying, but we have to duplicate what we have already
    // done and call it by a slightly different name if we want to have the same
    // function in multiple files
    double sum_over_t_edge_probability2 (
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
                double sum_term = sum_over_t_edge_probability2 (
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
}

using namespace Rcpp;
// for testing we will wrap and export this function so it is available in R.
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
