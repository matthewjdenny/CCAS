// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// ep
int ep(arma::vec intercepts, arma::mat coefficients, NumericVector latent_pos, int sender, int recipient, int actual_edge_value, arma::vec current_covariates, int interaction_pattern, bool using_coefficients);
RcppExport SEXP CCAS_ep(SEXP interceptsSEXP, SEXP coefficientsSEXP, SEXP latent_posSEXP, SEXP senderSEXP, SEXP recipientSEXP, SEXP actual_edge_valueSEXP, SEXP current_covariatesSEXP, SEXP interaction_patternSEXP, SEXP using_coefficientsSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::vec >::type intercepts(interceptsSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type coefficients(coefficientsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type latent_pos(latent_posSEXP);
    Rcpp::traits::input_parameter< int >::type sender(senderSEXP);
    Rcpp::traits::input_parameter< int >::type recipient(recipientSEXP);
    Rcpp::traits::input_parameter< int >::type actual_edge_value(actual_edge_valueSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type current_covariates(current_covariatesSEXP);
    Rcpp::traits::input_parameter< int >::type interaction_pattern(interaction_patternSEXP);
    Rcpp::traits::input_parameter< bool >::type using_coefficients(using_coefficientsSEXP);
    __result = Rcpp::wrap(ep(intercepts, coefficients, latent_pos, sender, recipient, actual_edge_value, current_covariates, interaction_pattern, using_coefficients));
    return __result;
END_RCPP
}
// lsms
int lsms(arma::vec unnormalized_discrete_distribution, int seed);
RcppExport SEXP CCAS_lsms(SEXP unnormalized_discrete_distributionSEXP, SEXP seedSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::vec >::type unnormalized_discrete_distribution(unnormalized_discrete_distributionSEXP);
    Rcpp::traits::input_parameter< int >::type seed(seedSEXP);
    __result = Rcpp::wrap(lsms(unnormalized_discrete_distribution, seed));
    return __result;
END_RCPP
}
