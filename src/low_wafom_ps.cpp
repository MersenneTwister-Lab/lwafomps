#include "low_wafom_ps.h"
#include <stdio.h>
#include <string.h>

namespace {
    struct lwps {
        int s;
        int d;
        double wafom;
        double discrepancy;
        uint32_t data[32 * 32];
    };
#include "lwps.h"
}

namespace LowWAFOMPointSet {
    //const double PointSet::TWO33 = 8589934592.0;
    const double PointSet::TWO53 = 9007199254740992.0;
    PointSet::PointSet()
    {
        data = NULL;
        set = NULL;
        start = 0;
        //end = 0;
        shift = 0;
        max = 0;
    }

    int PointSet::search(int s_value, int d_value)
    {
        if (data != NULL) {
            for (int i = 0; i < s; i++) {
                delete[] data[i];
            }
            delete[] data;
            data = NULL;
        }
        if (set != NULL) {
            delete[] set;
            set = NULL;
        }
        int last_index = -1;
        for (int i = 0; i < lwps_data_size; i++) {
            if (lwps_data[i].s == s_value) {
                if (lwps_data[i].d == d_value) {
                    last_index = i;
                    break;
                }
                if (lwps_data[i].d < d_value) {
                    last_index = i;
                }
            } else if (lwps_data[i].s > s_value) {
                break;
            }
        }
        //
        if (last_index >= 0) {
            s = lwps_data[last_index].s;
            d = lwps_data[last_index].d;
            data = new uint32_t *[s];
            for (int j = 0; j < s; j++) {
                data[j] = new uint32_t[d];
                for (int k = 0; k < d; k++) {
                    data[j][k] = lwps_data[last_index].data[k * s + j];
                }
            }
            reset();
            return d;
        }
        //
        return -1;
    }

    void PointSet::reset()
    {
        if (s == 0) {
            return;
        }
        if (set != NULL) {
            delete[] set;
        }
        set = new uint32_t[s];
        memset(set, 0, sizeof(uint32_t) * s);
        gray.reset();
        start = 0;
    }

    void PointSet::next()
    {
        gray.next();
        int bit = gray.index();
        for (int i = 0; i < s; i++) {
            set[i] ^= data[i][bit];
        }
    }
}
