// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// KendallNeighbour
NumericMatrix KendallNeighbour(NumericVector rank);
RcppExport SEXP rankdist_KendallNeighbour(SEXP rankSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type rank(rankSEXP);
    __result = Rcpp::wrap(KendallNeighbour(rank));
    return __result;
END_RCPP
}
// CayleyNeighbour
NumericMatrix CayleyNeighbour(NumericVector rank);
RcppExport SEXP rankdist_CayleyNeighbour(SEXP rankSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type rank(rankSEXP);
    __result = Rcpp::wrap(CayleyNeighbour(rank));
    return __result;
END_RCPP
}
// LogC
double LogC(NumericVector fai);
RcppExport SEXP rankdist_LogC(SEXP faiSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type fai(faiSEXP);
    __result = Rcpp::wrap(LogC(fai));
    return __result;
END_RCPP
}
// CWeightGivenPi
NumericVector CWeightGivenPi(NumericVector r1, NumericVector r2);
RcppExport SEXP rankdist_CWeightGivenPi(SEXP r1SEXP, SEXP r2SEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type r1(r1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type r2(r2SEXP);
    __result = Rcpp::wrap(CWeightGivenPi(r1, r2));
    return __result;
END_RCPP
}
// FindV
NumericMatrix FindV(NumericMatrix obs, NumericVector pi0);
RcppExport SEXP rankdist_FindV(SEXP obsSEXP, SEXP pi0SEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type obs(obsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type pi0(pi0SEXP);
    __result = Rcpp::wrap(FindV(obs, pi0));
    return __result;
END_RCPP
}
// LogC_Component
double LogC_Component(NumericVector fai);
RcppExport SEXP rankdist_LogC_Component(SEXP faiSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type fai(faiSEXP);
    __result = Rcpp::wrap(LogC_Component(fai));
    return __result;
END_RCPP
}
// Wtau
NumericMatrix Wtau(NumericMatrix obs, NumericVector pi0);
RcppExport SEXP rankdist_Wtau(SEXP obsSEXP, SEXP pi0SEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type obs(obsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type pi0(pi0SEXP);
    __result = Rcpp::wrap(Wtau(obs, pi0));
    return __result;
END_RCPP
}
// AllPerms
NumericMatrix AllPerms(int nobj);
RcppExport SEXP rankdist_AllPerms(SEXP nobjSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< int >::type nobj(nobjSEXP);
    __result = Rcpp::wrap(AllPerms(nobj));
    return __result;
END_RCPP
}
