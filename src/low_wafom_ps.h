#pragma once
#ifndef LOW_WAFORM_PS_H
#define LOW_WAFORM_PS_H
#include <stdint.h>
#include <stdio.h>
#include "gray.h"

namespace LowWAFOMPointSet {
    class PointSet {
    public:
        PointSet();
        ~PointSet() {
            if (data == NULL) {
                return;
            }
            for (int i = 0; i < s; i++) {
                delete[] data[i];
            }
            delete[] data;
        }
        int search(int s, int d);
        void reset();
        void next() ;
        /**
         *
         * @return r (0 < r < 1)
         */
        double operator[](int index) const {
            uint64_t x = set[index] ^ shift;
            x = x << (53 - 32);
            x |= 1;
            return static_cast<double>(x) / TWO53;
//            return (2.0 * static_cast<double>(set[index] ^ shift) + 1.0)
//                / TWO33;
        }
        bool hasNext() {
            return start < max;
        }
        void setShift(uint32_t value) {
            shift = value;
            reset();
        }
        int getS() {
            return s;
        }
        int getD() {
            return d;
        }
        double getWAFOM() {
            return wafom;
        }
        double getDiscrepancy() {
            return discrepancy;
        }
    private:
//        static const double TWO33;
        static const double TWO53;
        uint32_t **data;
        uint32_t *set;
        uint32_t shift;
        int start;
        int max;
        int s;
        int d;
        double wafom;
        double discrepancy;
        Gray gray;
    };
}
#endif // LOW_WAFORM_PS_H
