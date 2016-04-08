#include "mvnorm.h"
#include "kahan.h"
#include "inttypes.h"
#include <iostream>
#include <iomanip>
#include <math.h>
#include <fstream>
#include "erfinv.h"
#if !defined(NO_RCPP)
#include <Rcpp.h>
#endif
#include "low_wafom_ps.h"

using namespace std;
using namespace ErrorFunction;

namespace {
//const double Sqrt_PI = 1.7724538509055160272981674833411451827975494561224;
    const double Sqrt_PI = 1.772453850905516027298;
    /*
     * F_i(x) := \int_{-\infty }^{x} exp(-c_i t^2) dt is restated as follows
     * shc is precalculated shc = sqrt(0.5 * c)
     * (0.5 * Sqrt_PI / shc[i]) can be precalculated, but varsubF is called only in
     * precalcuation process.
     */
    void varsubF(int s, double result[], const double shc[], const double x[])
    {
        for (int i = 0; i < s; i++) {
            result[i] = (0.5 * Sqrt_PI / shc[i]) * (1.0 + erf(shc[i] * x[i]));
        }
    }
}

namespace LowWAFOMPointSet {
    PrecomputedData::PrecomputedData(int dimension_s)
    {
        s = dimension_s;
        upper = new double[s];
        calculatedUpper = new double[s];
        lower = new double[s];
        calculatedLower = new double[s];
        mean = new double[s];
        covarDiagC = new double[s];
        sqrtHalfC = new double[s];
        covarMatA_C = new double[s * s];
    }

    PrecomputedData::~PrecomputedData()
    {
        if (upper != NULL) {
            delete[] upper;
        }
        if (calculatedUpper != NULL) {
            delete[] calculatedUpper;
        }
        if (lower != NULL) {
            delete[] lower;
        }
        if (calculatedLower != NULL) {
            delete[] calculatedLower;
        }
        if (mean != NULL) {
            delete[] mean;
        }
        if (covarDiagC != NULL) {
            delete[] covarDiagC;
        }
        if (sqrtHalfC != NULL) {
            delete[] sqrtHalfC;
        }
        if (covarMatA_C != NULL) {
            delete[] covarMatA_C;
        }
    }

    void PrecomputedData::print(ostream& os)
    {
        os.precision(20);
        os << "s:" << dec << s << endl;
        os << "productFactor:" << productFactor << endl;
        for (int i = 0; i < s; i++) {
            os << "upper[" << dec << i << "]:" << upper[i] << endl;
        }
        for (int i = 0; i < s; i++) {
            os << "lower[" << dec << i << "]:" << lower[i] << endl;
        }
        for (int i = 0; i < s; i++) {
            os << "mean[" << dec << i << "]:" << mean[i] << endl;
        }
        for (int i = 0; i < s; i++) {
            os << "covarDiagC[" << dec << i << "]:" << covarDiagC[i] << endl;
        }
        for (int i = 0; i < s; i++) {
            os << "sqrtHalfC[" << dec << i << "]:" << sqrtHalfC[i] << endl;
        }
        for (int i = 0; i < s; i++) {
            os << "calculatedUpper[" << dec << i << "]:"
               << calculatedUpper[i] << endl;
        }
        for (int i = 0; i < s; i++) {
            os << "calculatedLower[" << dec << i << "]:"
               << calculatedLower[i] << endl;
        }
        for (int i = 0; i < s; i++) {
            for (int j = 0; j < s; j++) {
                os << "covarMatA_C[" << i << "," << j << "]:"
                   << covarMatA_C[i * s + j] << endl;
            }
        }
    }

    /**
     * PreComputeation
     * This will be usefull for randomized MCQMC
     */
    void MVNorm_LWAFOMPS::precompute(PrecomputedData& result,
                                     int s_dimension,
                                     const double lower[],
                                     const double upper[],
                                     const double mean[],
                                     const double coval[])
    {
        // result.s = s_dimension;   // already set.
        int s = s_dimension;
        for (int i = 0; i < s; i++) {
            result.upper[i] = upper[i] - mean[i];
            result.lower[i] = lower[i] - mean[i];
            result.mean[i] = 0;
        }
        // Set Up result.covarDiagC
        Kahan sum;
        for (int i = 0; i < s; i++) {
            sum.clear();
            for (int  j = 0; j < s; j++) {
                if (j != i) {
                    sum.add(coval[i * s + j]);
                }
            }
            result.covarDiagC[i] = coval[i * s + i] - sum.get();
        }
        // Set Up result.covarMatA_C
        for (int i = 0; i < s; ++i) {
            for (int j = 0; j < s; ++j) {
                if ( i != j ) {
                    result.covarMatA_C[i * s + j] = coval[i * s + j];
                } else {
                    result.covarMatA_C[i * s + j]
                        = coval[i * s + j] - result.covarDiagC[i];
                }
            }
        }
        // Set Up result.sqrtHalf
        for (int i = 0; i < s; i++) {
            result.sqrtHalfC[i] = sqrt(0.5 * result.covarDiagC[i]);
        }
        // Set Up result.productFactor
        double u[s];
        double l[s];
        varsubF(s, u, result.sqrtHalfC, result.upper);
        varsubF(s, l, result.sqrtHalfC, result.lower);
        result.productFactor = 1.0;
        for (int i = 0; i < s; i++) {
            if (result.lower[i] == -INFINITY) {
                result.productFactor *= u[i];
            } else {
                result.productFactor *= u[i] - l[i];
            }
        }
        // Set Up result.calculatedUpper
        for (int i = 0; i < s; i++) {
            result.calculatedUpper[i]
                = erf(result.sqrtHalfC[i] * result.upper[i]);
            result.calculatedLower[i]
                = erf(result.sqrtHalfC[i] * result.lower[i]);
        }
    }

    MVNorm_LWAFOMPS::MVNorm_LWAFOMPS(const PrecomputedData& precomputedData)
    {
        data = new PrecomputedData(precomputedData.s);
        data->s = precomputedData.s;
        data->productFactor = precomputedData.productFactor;
        for (int i = 0; i < data->s; i++) {
            data->upper[i] = precomputedData.upper[i];
            data->lower[i] = precomputedData.lower[i];
            data->mean[i] = precomputedData.mean[i];
            data->covarDiagC[i] = precomputedData.covarDiagC[i];
            data->sqrtHalfC[i] = precomputedData.sqrtHalfC[i];
            data->calculatedUpper[i] = precomputedData.calculatedUpper[i];
            data->calculatedLower[i] = precomputedData.calculatedLower[i];
        }
        for (int i = 0; i < data->s; i++) {
            for (int j = 0; j < data->s; j++) {
                data->covarMatA_C[i * data->s + j]
                    = precomputedData.covarMatA_C[i * data->s + j];
            }
        }
    }

    MVNorm_LWAFOMPS::~MVNorm_LWAFOMPS()
    {
        if (data != NULL) {
            delete data;
        }
    }
#if 0
    /**
     *  E_i(F_i(b_i)*x_i) is expressed as follows
     * @param result outputs
     * @param shc sqrt(c * 0.5)
     * @param esb erf(shc * b)
     * @param x points
     */
    void  MVNorm_LWAFOMPS::ArgOfExp(double result[], PointSet& ps)
    {
        int s = data->s;
        for (int i = 0; i < s; i++) {
            result[i] =  (1.0 / data->sqrtHalfC[i])
                * erfinv((1.0 + data->calculatedUpper[i]) * ps[i] - 1.0);
        }
    }
#endif
    /**
     * integrand function
     * multivarite normal distribution
     * @param shc sqrt(c * 0.5)
     * @param esb erf(shc * b)
     */
    double MVNorm_LWAFOMPS::integrand(PointSet& ps)
    {
        int s = data->s;
        Kahan sum;
        double memo[s];
        //ArgOfExp(memo,  ps);
        //
        for (int i = 0; i < s; i++) {
            if (data->lower[i] == -INFINITY) {
                memo[i] = (1.0 / data->sqrtHalfC[i])
                    * erfinv((1.0 + data->calculatedUpper[i]) * ps[i] - 1.0);
            } else {
                memo[i] = (1.0 / data->sqrtHalfC[i])
                    * erfinv(
                        (data->calculatedUpper[i]
                         - data->calculatedLower[i]) * ps[i]
                        + data->calculatedLower[i]);
            }
        }
        //
        for (int i = 0; i < s; i++) {
            sum.add(memo[i] * data->covarMatA_C[i * s + i] * memo[i]);
        }

        for (int i = 0; i < s; i++) {
            for (int j = i + 1; j < s; j++) {
                sum.add(2.0 * memo[i] * data->covarMatA_C[i * s + j] * memo[j]);
            }
        }
        return exp(-0.5 * sum.get());
    }

    double MVNorm_LWAFOMPS::integrate(int d, PointSet& ps, uint32_t shift)
    {
#if !defined(NO_RCPP)
        Rcpp::checkUserInterrupt();
#endif
        Kahan QMCvalue;
        ps.setShift(shift);
        int32_t max = INT32_C(1) << d;
        for (int a = 0; a < max; ++a) {
#if !defined(NO_RCPP)
            if (a % 10000 == 0) {
                Rcpp::checkUserInterrupt();
            }
#endif
            QMCvalue.add(integrand(ps));
            ps.next();
        }
        double N = static_cast<double>(max);
        return data->productFactor * QMCvalue.get() / N;
    }
}
