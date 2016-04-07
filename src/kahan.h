#ifndef KAHAN_HPP
#define KAHAN_HPP

/*
 * Kahan summation algorithm
 */
class Kahan {
public:
    Kahan() {
        sum = 0.0;
        c = 0.0;
    }
    void clear() {
        sum = 0.0;
        c = 0.0;
    }
    void add(double x) {
        double y = x - c;
        double t = sum + y;
        c = (t - sum) - y;
        sum = t;
    }
    double get() {
        return sum;
    }
private:
    double c;
    double sum;
};
#endif // KAHAN_HPP
