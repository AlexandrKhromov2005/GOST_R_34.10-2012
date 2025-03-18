#ifndef POINT_HPP
#define POINT_HPP

#include <gmpxx.h>
#include <stdexcept>
#include "algorithms_for_primes.hpp"
#include <vector>
#include <map>
#include <cmath>
#include <iostream>

class Curve;

class Point {
private:
    mpz_class x;
    mpz_class y;
    bool is_infinity;
    Curve* crv;
    mpz_class order;

public:
    Point(Curve* crv);
    Point(const mpz_class& x, const mpz_class& y, Curve* crv, bool is_infinity);
    Point(const Point& other) = default; 
    
    mpz_class get_x() const;
    mpz_class get_y() const;
    bool isInfinity() const;
    Curve* get_curve() const;

    Point operator+(const Point& other) const;
    Point operator*(const mpz_class &k) const;
    bool operator==(const Point& other) const;
    mpz_class refine_order(mpz_class candidate) const;
    void calculate_order();
    std::string to_string() const;
    mpz_class find_order_bruteforce(const mpz_class& N) const;
    void set_order(mpz_class ord);
    mpz_class get_order() const;
};

#endif // POINT_HPP