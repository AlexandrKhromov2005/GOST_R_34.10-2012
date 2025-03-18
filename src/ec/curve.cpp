#include "curve.hpp"
#include <iostream>
#include <set>

Curve::Curve(mpz_class a, mpz_class b, mpz_class p) : a(a), b(b), p(p) {}
Curve::Curve() : a(0), b(0), p(0) {}

mpz_class Curve::get_a() { return a; }
mpz_class Curve::get_b() { return b; }
mpz_class Curve::get_p() { return p; }

void Curve::find_points() {
    points.push_back(Point(this)); 

    for (mpz_class x = 0; x < p; ++x) {
        mpz_class y_pow_2 = (x * x * x + a * x + b) % p;
        if (y_pow_2 == 0) {
            Point pt(x, 0, this, false);
            points.push_back(pt);
            continue;
        }

        if (legendre_symbol(y_pow_2, p) == 1) {
            Point pt1(x, mod_sqrt(y_pow_2, p), this, false);
            points.push_back(pt1);

            Point pt2(x, p - mod_sqrt(y_pow_2, p), this, false);
            points.push_back(pt2);
        }
    }
}

void Curve::print_points() {
    mpz_class i = 0;
    for (const Point& point : points) {
        if (point.isInfinity()) {
            std::cout << "Infinity point" << std::endl;
            ++i;
            continue;
        }
        std::cout << "Point " << i << ": (" << point.get_x() << " , " << point.get_y() << ")"
                  << " Order: " << point.get_order() << std::endl;
        ++i;
    }

    printPrimeSubgroups();
}

bool Curve::is_on_curve(const Point& P) const {
    if (P.isInfinity()) return true; 
    mpz_class x = P.get_x();
    mpz_class y = P.get_y();

    mpz_class lhs = (y * y) % p;         
    mpz_class rhs = (x * x * x + a * x + b) % p; 
    return lhs == rhs;
}


mpz_class Curve::get_group_order() {
    return mpz_class(static_cast<unsigned long>(points.size())); 
}

void Curve::printPrimeSubgroups() {
    mpz_class order = get_group_order();
    std::vector<mpz_class> prime_factors = factorize(order); 
    std::set<mpz_class> unique_factors(prime_factors.begin(), prime_factors.end());
    std::vector<mpz_class> factors(unique_factors.begin(), unique_factors.end());

    std::set<std::string> processed_subgroups;

    for (const mpz_class& factor : factors) {
        for (Point &point : points) {
            if (point.get_order() == factor) {
                std::vector<Point> subgroup;
                for (mpz_class k = 0; k < factor; ++k) {
                    subgroup.push_back(point * k);
                }

                std::vector<Point> sorted_subgroup = subgroup;
                std::sort(sorted_subgroup.begin(), sorted_subgroup.end(),
                         [](const Point& a, const Point& b) {
                             if (a.get_x() != b.get_x()) {
                                 return a.get_x() < b.get_x();
                             }
                             return a.get_y() < b.get_y();
                         });

                std::string key;
                for (const Point& p : sorted_subgroup) {
                    key += "(" + p.get_x().get_str() + "," + p.get_y().get_str() + ")";
                }

                if (!processed_subgroups.count(key)) {
                    processed_subgroups.insert(key);

                    std::cout << "Subgroup of prime order " << factor << ":\n";
                    for (const Point& p : subgroup) {
                        if (!p.isInfinity()){
                            std::cout << "Point: (" << p.get_x() << ", " << p.get_y() << ")\n";
                        } else {
                            std::cout << "infinite point" << std::endl;
                        }
                    }
                    std::cout << std::endl;
                }
            }
        }
    }
}


Point Curve::find_point_of_order(const mpz_class &q) {
    mpz_class group_order = get_group_order();
    
    // Проверка: q должен делить порядок группы
    if (group_order % q != 0) {
        throw std::runtime_error("q не делит порядок группы точек");
    }

    // Проходим по всем точкам, выбираем случайную (не бесконечную)
    for (auto &R : points) {
        if (R.isInfinity()) {
            continue;
        }
        // Вычисляем P = (group_order / q) * R
        Point P = R * (group_order / q);
        
        // Если P не равна бесконечной точке и q * P = бесконечность, то нашли нужную точку
        if (!P.isInfinity() && (P * q).isInfinity()) {
            return P;
        }
    }
    throw std::runtime_error("Точка с требуемым порядком не найдена");
}