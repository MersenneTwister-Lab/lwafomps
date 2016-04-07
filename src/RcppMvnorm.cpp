#include <Rcpp.h>
#include <math.h>
#include "mvnorm.h"
#include "low_wafom_ps.h"
#include "kahan.h"

using namespace std;
using namespace Rcpp;
using namespace LowWAFOMPointSet;

//Pre-computation for Multivariate Normal Distributuion
//@param lower vector
//@param upper vector
//@param mean vector
//@param coval Variance-covariance matrix
//@return pre-computed data
// [[Rcpp::export(rng = false)]]
List rcppPrecompute(NumericVector lower,
                    NumericVector upper,
                    NumericVector mean,
                    NumericMatrix coval)
{
    //Rprintf("rcppPrecompute start\n");
    int s = lower.length();
    if (s != upper.length() || s != mean.length()
        || s != coval.nrow() || s != coval.ncol()) {
        Rcpp::stop("dimension mismatch!");
    }
    // check if coval is symmetric matrix
    for (int i = 0; i < coval.nrow(); i++) {
        for (int j = i + 1; j < coval.ncol(); j++) {
            if (coval(i, j) != coval(j,i)) {
                Rcpp::stop("coval is not symmetric");
            }
        }
    }
    // check if coval is diagonally dominant matrix
    Kahan sum;
    for (int i = 0; i < coval.nrow(); i++) {
        sum.clear();
        for (int j = 0; j < coval.ncol(); j++) {
            if (i == j) {
                continue;
            }
            sum.add(fabs(coval(i, j)));
        }
        if (fabs(coval(i, i)) < sum.get()) {
            Rcpp::stop("coval is not diagonally dominant");
        }
    }
    double cppLower[s];
    double cppUpper[s];
    double cppMean[s];
    for (int i = 0; i < s; i++) {
        cppLower[i] = lower[i];
        cppUpper[i] = upper[i];
        cppMean[i] = mean[i];
    }
    double cppCoval[s * s];
    for (int i = 0; i < s; i++) {
        for (int j = 0; j < s; j++) {
            cppCoval[i * s + j] = coval(i, j);
        }
    }
    PrecomputedData cppData(s);
    MVNorm_LWAFOMPS::precompute(cppData, s, cppLower,
                                cppUpper, cppMean, cppCoval);
    NumericVector covarDiagC(s);
    NumericVector sqrtHalfC(s);
    NumericVector calculatedUpper(s);
    NumericMatrix covarMatAC(s,s);
    for (int i = 0; i < s; i++) {
        lower[i] = cppData.lower[i];
        upper[i] = cppData.upper[i];
        covarDiagC[i] = cppData.covarDiagC[i];
        sqrtHalfC[i] = cppData.sqrtHalfC[i];
        calculatedUpper[i] = cppData.calculatedUpper[i];
    }
    for (int i = 0; i < s; i++) {
        for (int j = 0; j < s; j++) {
            covarMatAC(i, j) = cppData.covarMatA_C[i * s + j];
        }
    }
    List data = List::create(Named("dimension")=cppData.s,
                             Named("productFactor")=cppData.productFactor,
                             Named("lower")=lower,
                             Named("upper")=upper,
                             Named("covarDiag")=covarDiagC,
                             Named("sqrtHalfC")=sqrtHalfC,
                             Named("calculatedUpper")=calculatedUpper,
                             Named("covarMatAC")=covarMatAC);
    return data;
}

//Multivariate Normal Distributuion
//@param rcpData pre-computed data
//@param shift random shift
//@return integrated value
// [[Rcpp::export(rng = false)]]
double rcppMvnorm(List rcpData, int shift)
{
    //Rprintf("rcppMvnorm start\n");
    int s = as<int>(rcpData[0]);
    double productFactor = as<double>(rcpData[1]);
    NumericVector lower = as<NumericVector>(rcpData[2]);
    NumericVector upper = as<NumericVector>(rcpData[3]);
    NumericVector covarDiagC = as<NumericVector>(rcpData[4]);
    NumericVector sqrtHalfC = as<NumericVector>(rcpData[5]);
    NumericVector calculatedUpper = as<NumericVector>(rcpData[6]);
    NumericMatrix covarMatAC = as<NumericMatrix>(rcpData[7]);
    if (s != lower.length()
        || s != upper.length()
        || s != covarDiagC.length()
        || s != sqrtHalfC.length()
        || s != calculatedUpper.length()
        || s != covarMatAC.nrow()
        || s != covarMatAC.ncol()) {
        Rcpp::stop("dimension mismatch!");
    }

    PrecomputedData cppData(s);
    cppData.productFactor = productFactor;
    for (int i = 0; i < s; i++) {
        cppData.lower[i] = lower[i];
        cppData.upper[i] = upper[i];
        cppData.covarDiagC[i] = covarDiagC[i];
        cppData.sqrtHalfC[i] = sqrtHalfC[i];
        cppData.calculatedUpper[i] = calculatedUpper[i];
    }
    for (int i = 0; i < s; i++) {
        for (int j = 0; j < s; j++) {
            cppData.covarMatA_C[i * s + j] = covarMatAC(i,j);
        }
    }
    MVNorm_LWAFOMPS mvnorm(cppData);
    //int d = 20;
    int d = 1000;
    PointSet ps;
    int found_d = ps.search(s, d);
    if (found_d < 0) {
        Rcpp::stop("can't find Low WAFOM Point Set");
    }
    double QMCvalue = mvnorm.integrate(found_d, ps, shift);
    //Rprintf("QMCvalue = %f\n", QMCvalue);
    return QMCvalue;
}
