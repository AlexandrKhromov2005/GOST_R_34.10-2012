#ifndef CURVE_HPP
#define CURVE_HPP

#include <vector>
#include "point.hpp"
#include "algorithms_for_primes.hpp"

class Curve {
private:
    mpz_class a, b, p;

public:
    std::vector<Point> points; 
    std::vector<mpz_class> orders;

    Curve(mpz_class a, mpz_class b, mpz_class p);
    Curve();

    mpz_class get_a();
    mpz_class get_b();
    mpz_class get_p();

    Point find_point_of_order(const mpz_class &q);

    bool is_on_curve(const Point& P) const;

    void find_points();
    void print_points();
    mpz_class get_group_order();
    void printPrimeSubgroups();
};

#endif
