#ifndef GRAY_HPP
#define GRAY_HPP
#include <stdint.h>

class Gray {
public:
    Gray() {
        count = 0;
        pre = 0;
        gray = 0;
    }
    uint32_t next() {
        count++;
        pre = gray;
        gray = count ^ (count >> 1);
        return gray;
    }
    void reset() {
        count = 0;
        pre = 0;
        gray = 0;
    }
    int index() {
        return bitPos(pre ^ gray);
    }
private:
    uint32_t count;
    uint32_t gray;
    uint32_t pre;
    int ones32(uint32_t x) {
        x -= ((x >> 1) & 0x55555555);
        x = (((x >> 2) & 0x33333333) + (x & 0x33333333));
        x = (((x >> 4) + x) & 0x0f0f0f0f);
        x += (x >> 8);
        x += (x >> 16);
        return(x & 0x0000003f);
    }
    int bitPos(uint32_t x) {
        x |= (x >> 1);
        x |= (x >> 2);
        x |= (x >> 4);
        x |= (x >> 8);
        x |= (x >> 16);
        return(ones32(x >> 1));
    }
};
#endif
