#pragma once
#ifndef MVNORM_H
#define MVNORM_H
#include <iostream>
#include "low_wafom_ps.h"

namespace LowWAFOMPointSet {
    class PrecomputedData {
    public:
        PrecomputedData(int s);
        ~PrecomputedData();
        void print(std::ostream& os);
        int s;
        double productFactor;
        double *upper;
        double *lower;
        double *mean;
        double *covarDiagC;
        double *sqrtHalfC;
        double *calculatedUpper;
        double *calculatedLower;
        double *covarMatA_C;
    };

    class MVNorm_LWAFOMPS {
    public:
        static void precompute(PrecomputedData& result,
                               int s_dimension,
                               const double lower[],
                               const double upper[],
                               const double mean[],
                               const double coval[]);
        MVNorm_LWAFOMPS(const PrecomputedData& precomputedData);
        ~MVNorm_LWAFOMPS();
        /*
         *
         */
        double integrate(int d, PointSet& ps, uint32_t shift = 0);
    private:
        //void varsubE(double result[], const double t[]);
        //void  ArgOfExp(double result[], PointSet& ps);
        double integrand(PointSet& ps);
        PrecomputedData *data;
    };
}
#endif
