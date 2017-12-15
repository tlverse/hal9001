// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "hal9001_types.h"
#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// index_first_copy
IntegerVector index_first_copy(const MSpMat& X);
RcppExport SEXP _hal9001_index_first_copy(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const MSpMat& >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(index_first_copy(X));
    return rcpp_result_gen;
END_RCPP
}
// apply_copy_map
SpMat apply_copy_map(const MSpMat& X, const List& copy_map);
RcppExport SEXP _hal9001_apply_copy_map(SEXP XSEXP, SEXP copy_mapSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const MSpMat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const List& >::type copy_map(copy_mapSEXP);
    rcpp_result_gen = Rcpp::wrap(apply_copy_map(X, copy_map));
    return rcpp_result_gen;
END_RCPP
}
// lassi_predict
NumericVector lassi_predict(const MSpMat X, const NumericVector beta, double intercept);
RcppExport SEXP _hal9001_lassi_predict(SEXP XSEXP, SEXP betaSEXP, SEXP interceptSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const MSpMat >::type X(XSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< double >::type intercept(interceptSEXP);
    rcpp_result_gen = Rcpp::wrap(lassi_predict(X, beta, intercept));
    return rcpp_result_gen;
END_RCPP
}
// soft_threshold
double soft_threshold(double beta, double lambda);
RcppExport SEXP _hal9001_soft_threshold(SEXP betaSEXP, SEXP lambdaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    rcpp_result_gen = Rcpp::wrap(soft_threshold(beta, lambda));
    return rcpp_result_gen;
END_RCPP
}
// X_t_resid
double X_t_resid(const MSpMat& X, const NumericVector& resids, int j, double xscale_j, double xcenter_j, double resid_sum);
RcppExport SEXP _hal9001_X_t_resid(SEXP XSEXP, SEXP residsSEXP, SEXP jSEXP, SEXP xscale_jSEXP, SEXP xcenter_jSEXP, SEXP resid_sumSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const MSpMat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type resids(residsSEXP);
    Rcpp::traits::input_parameter< int >::type j(jSEXP);
    Rcpp::traits::input_parameter< double >::type xscale_j(xscale_jSEXP);
    Rcpp::traits::input_parameter< double >::type xcenter_j(xcenter_jSEXP);
    Rcpp::traits::input_parameter< double >::type resid_sum(resid_sumSEXP);
    rcpp_result_gen = Rcpp::wrap(X_t_resid(X, resids, j, xscale_j, xcenter_j, resid_sum));
    return rcpp_result_gen;
END_RCPP
}
// get_new_beta
double get_new_beta(const MSpMat& X, const NumericVector& resids, int j, double beta_j, double xscale_j, double xcenter_j, double resid_sum);
RcppExport SEXP _hal9001_get_new_beta(SEXP XSEXP, SEXP residsSEXP, SEXP jSEXP, SEXP beta_jSEXP, SEXP xscale_jSEXP, SEXP xcenter_jSEXP, SEXP resid_sumSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const MSpMat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type resids(residsSEXP);
    Rcpp::traits::input_parameter< int >::type j(jSEXP);
    Rcpp::traits::input_parameter< double >::type beta_j(beta_jSEXP);
    Rcpp::traits::input_parameter< double >::type xscale_j(xscale_jSEXP);
    Rcpp::traits::input_parameter< double >::type xcenter_j(xcenter_jSEXP);
    Rcpp::traits::input_parameter< double >::type resid_sum(resid_sumSEXP);
    rcpp_result_gen = Rcpp::wrap(get_new_beta(X, resids, j, beta_j, xscale_j, xcenter_j, resid_sum));
    return rcpp_result_gen;
END_RCPP
}
// find_lambda_max
double find_lambda_max(const MSpMat& X, const NumericVector& y, const NumericVector& xscale, const NumericVector& xcenter);
RcppExport SEXP _hal9001_find_lambda_max(SEXP XSEXP, SEXP ySEXP, SEXP xscaleSEXP, SEXP xcenterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const MSpMat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type xscale(xscaleSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type xcenter(xcenterSEXP);
    rcpp_result_gen = Rcpp::wrap(find_lambda_max(X, y, xscale, xcenter));
    return rcpp_result_gen;
END_RCPP
}
// equal_double
bool equal_double(double x, double y);
RcppExport SEXP _hal9001_equal_double(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(equal_double(x, y));
    return rcpp_result_gen;
END_RCPP
}
// update_coord
double update_coord(const MSpMat& X, NumericVector& resids, NumericVector& beta, double lambda, int j, const NumericVector& xscale, const NumericVector& xcenter, double& resid_sum, bool center);
RcppExport SEXP _hal9001_update_coord(SEXP XSEXP, SEXP residsSEXP, SEXP betaSEXP, SEXP lambdaSEXP, SEXP jSEXP, SEXP xscaleSEXP, SEXP xcenterSEXP, SEXP resid_sumSEXP, SEXP centerSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const MSpMat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type resids(residsSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< int >::type j(jSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type xscale(xscaleSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type xcenter(xcenterSEXP);
    Rcpp::traits::input_parameter< double& >::type resid_sum(resid_sumSEXP);
    Rcpp::traits::input_parameter< bool >::type center(centerSEXP);
    rcpp_result_gen = Rcpp::wrap(update_coord(X, resids, beta, lambda, j, xscale, xcenter, resid_sum, center));
    return rcpp_result_gen;
END_RCPP
}
// update_coords
int update_coords(const MSpMat& X, NumericVector& resids, NumericVector& beta, double lambda, const NumericVector& xscale, const NumericVector& xcenter, NumericVector& intercept, bool active_set, bool center);
RcppExport SEXP _hal9001_update_coords(SEXP XSEXP, SEXP residsSEXP, SEXP betaSEXP, SEXP lambdaSEXP, SEXP xscaleSEXP, SEXP xcenterSEXP, SEXP interceptSEXP, SEXP active_setSEXP, SEXP centerSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const MSpMat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type resids(residsSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type xscale(xscaleSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type xcenter(xcenterSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type intercept(interceptSEXP);
    Rcpp::traits::input_parameter< bool >::type active_set(active_setSEXP);
    Rcpp::traits::input_parameter< bool >::type center(centerSEXP);
    rcpp_result_gen = Rcpp::wrap(update_coords(X, resids, beta, lambda, xscale, xcenter, intercept, active_set, center));
    return rcpp_result_gen;
END_RCPP
}
// lassi_fit_cd
int lassi_fit_cd(const MSpMat& X, NumericVector& resids, NumericVector& beta, double lambda, int nsteps, const NumericVector& xscale, const NumericVector& xcenter, NumericVector& intercept, bool active_set, bool center);
RcppExport SEXP _hal9001_lassi_fit_cd(SEXP XSEXP, SEXP residsSEXP, SEXP betaSEXP, SEXP lambdaSEXP, SEXP nstepsSEXP, SEXP xscaleSEXP, SEXP xcenterSEXP, SEXP interceptSEXP, SEXP active_setSEXP, SEXP centerSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const MSpMat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type resids(residsSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< int >::type nsteps(nstepsSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type xscale(xscaleSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type xcenter(xcenterSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type intercept(interceptSEXP);
    Rcpp::traits::input_parameter< bool >::type active_set(active_setSEXP);
    Rcpp::traits::input_parameter< bool >::type center(centerSEXP);
    rcpp_result_gen = Rcpp::wrap(lassi_fit_cd(X, resids, beta, lambda, nsteps, xscale, xcenter, intercept, active_set, center));
    return rcpp_result_gen;
END_RCPP
}
// make_basis_list
List make_basis_list(const NumericMatrix& X_sub, const NumericVector& cols);
RcppExport SEXP _hal9001_make_basis_list(SEXP X_subSEXP, SEXP colsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type X_sub(X_subSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type cols(colsSEXP);
    rcpp_result_gen = Rcpp::wrap(make_basis_list(X_sub, cols));
    return rcpp_result_gen;
END_RCPP
}
// meets_basis
bool meets_basis(const NumericMatrix& X, const int row_num, const IntegerVector& cols, const NumericVector& cutoffs);
RcppExport SEXP _hal9001_meets_basis(SEXP XSEXP, SEXP row_numSEXP, SEXP colsSEXP, SEXP cutoffsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const int >::type row_num(row_numSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type cols(colsSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type cutoffs(cutoffsSEXP);
    rcpp_result_gen = Rcpp::wrap(meets_basis(X, row_num, cols, cutoffs));
    return rcpp_result_gen;
END_RCPP
}
// evaluate_basis
void evaluate_basis(const List& basis, const NumericMatrix& X, SpMat& x_basis, int basis_col);
RcppExport SEXP _hal9001_evaluate_basis(SEXP basisSEXP, SEXP XSEXP, SEXP x_basisSEXP, SEXP basis_colSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const List& >::type basis(basisSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type X(XSEXP);
    Rcpp::traits::input_parameter< SpMat& >::type x_basis(x_basisSEXP);
    Rcpp::traits::input_parameter< int >::type basis_col(basis_colSEXP);
    evaluate_basis(basis, X, x_basis, basis_col);
    return R_NilValue;
END_RCPP
}
// make_design_matrix
SpMat make_design_matrix(const NumericMatrix& X, const List& blist);
RcppExport SEXP _hal9001_make_design_matrix(SEXP XSEXP, SEXP blistSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const List& >::type blist(blistSEXP);
    rcpp_result_gen = Rcpp::wrap(make_design_matrix(X, blist));
    return rcpp_result_gen;
END_RCPP
}
// asdgCMatrix_
SEXP asdgCMatrix_(SEXP XX_);
RcppExport SEXP _hal9001_asdgCMatrix_(SEXP XX_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type XX_(XX_SEXP);
    rcpp_result_gen = Rcpp::wrap(asdgCMatrix_(XX_));
    return rcpp_result_gen;
END_RCPP
}
// non_zeros
IntegerVector non_zeros(const MSpMat& X);
RcppExport SEXP _hal9001_non_zeros(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const MSpMat& >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(non_zeros(X));
    return rcpp_result_gen;
END_RCPP
}
// get_pnz
NumericVector get_pnz(const MSpMat& X);
RcppExport SEXP _hal9001_get_pnz(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const MSpMat& >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(get_pnz(X));
    return rcpp_result_gen;
END_RCPP
}
// get_xscale
NumericVector get_xscale(const MSpMat& X, const NumericVector& xcenter);
RcppExport SEXP _hal9001_get_xscale(SEXP XSEXP, SEXP xcenterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const MSpMat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type xcenter(xcenterSEXP);
    rcpp_result_gen = Rcpp::wrap(get_xscale(X, xcenter));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_hal9001_index_first_copy", (DL_FUNC) &_hal9001_index_first_copy, 1},
    {"_hal9001_apply_copy_map", (DL_FUNC) &_hal9001_apply_copy_map, 2},
    {"_hal9001_lassi_predict", (DL_FUNC) &_hal9001_lassi_predict, 3},
    {"_hal9001_soft_threshold", (DL_FUNC) &_hal9001_soft_threshold, 2},
    {"_hal9001_X_t_resid", (DL_FUNC) &_hal9001_X_t_resid, 6},
    {"_hal9001_get_new_beta", (DL_FUNC) &_hal9001_get_new_beta, 7},
    {"_hal9001_find_lambda_max", (DL_FUNC) &_hal9001_find_lambda_max, 4},
    {"_hal9001_equal_double", (DL_FUNC) &_hal9001_equal_double, 2},
    {"_hal9001_update_coord", (DL_FUNC) &_hal9001_update_coord, 9},
    {"_hal9001_update_coords", (DL_FUNC) &_hal9001_update_coords, 9},
    {"_hal9001_lassi_fit_cd", (DL_FUNC) &_hal9001_lassi_fit_cd, 10},
    {"_hal9001_make_basis_list", (DL_FUNC) &_hal9001_make_basis_list, 2},
    {"_hal9001_meets_basis", (DL_FUNC) &_hal9001_meets_basis, 4},
    {"_hal9001_evaluate_basis", (DL_FUNC) &_hal9001_evaluate_basis, 4},
    {"_hal9001_make_design_matrix", (DL_FUNC) &_hal9001_make_design_matrix, 2},
    {"_hal9001_asdgCMatrix_", (DL_FUNC) &_hal9001_asdgCMatrix_, 1},
    {"_hal9001_non_zeros", (DL_FUNC) &_hal9001_non_zeros, 1},
    {"_hal9001_get_pnz", (DL_FUNC) &_hal9001_get_pnz, 1},
    {"_hal9001_get_xscale", (DL_FUNC) &_hal9001_get_xscale, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_hal9001(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
