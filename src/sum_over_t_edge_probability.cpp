// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(BH)]]

#include <RcppArmadillo.h>

using std::log;
using std::exp;
using std::max;
using std::abs;

namespace mjd {
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
}

using namespace Rcpp;
// for testing we will wrap and export this function so it is available in R.
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
