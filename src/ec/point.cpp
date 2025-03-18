#include "point.hpp"
#include "curve.hpp"
#include "algorithms_for_primes.hpp"
#include <set>

// Конструкторы
Point::Point(Curve* crv) : x(0), y(0), is_infinity(true), crv(crv) {}
Point::Point(const mpz_class& x, const mpz_class& y, Curve* crv, bool is_infinity) 
    : x(mod(x, crv->get_p())), 
      y(mod(y, crv->get_p())), 
      is_infinity(false), 
      crv(crv) {}

Point::Point(const Point& other)
    : x(other.x),
      y(other.y),
      is_infinity(other.is_infinity),
      crv(other.crv),
      order(other.order) {}

mpz_class Point::get_x() const { return x; }
mpz_class Point::get_y() const { return y; }
bool Point::isInfinity() const { return is_infinity; }
Curve* Point::get_curve() const { return crv; }

bool Point::operator==(const Point& other) const {
    if (isInfinity() && other.isInfinity()) return true;
    if (isInfinity() || other.isInfinity()) return false;
    return (x == other.x) && (y == other.y) && (crv == other.crv);
}

Point Point::operator+(const Point& other) const {
    if (crv != other.crv)
        throw std::invalid_argument("Points are on different curves");

    // Если одна из точек – нейтральный элемент, возвращаем другую
    if (this->isInfinity()) return other;
    if (other.isInfinity()) return *this;

    mpz_class p = crv->get_p();
    mpz_class x1 = x, y1 = y;
    mpz_class x2 = other.x, y2 = other.y;

    if (x1 == x2 && mod(y1 + y2, p) == 0) {
        return Point(crv); 
    }

    mpz_class lambda;
    if (*this == other) { 
        if (y1 == 0) {
            return Point(crv);
        }
        mpz_class numerator = mod(3 * mod(x1 * x1, p) + crv->get_a(), p);
        mpz_class denominator = mod(2 * y1, p);
        mpz_class inv_denominator = mod_inverse(denominator, p);
        lambda = mod(numerator * inv_denominator, p);
    } else { 
        mpz_class numerator = mod(y2 - y1, p);
        mpz_class denominator = mod(x2 - x1, p);
        mpz_class inv_denominator = mod_inverse(denominator, p);
        lambda = mod(numerator * inv_denominator, p);
    }

    mpz_class x3 = mod(lambda * lambda - x1 - x2, p);
    mpz_class y3 = mod(lambda * (x1 - x3) - y1, p);


    return Point(x3, y3, crv, false);
}



Point Point::operator*(const mpz_class &k) const {
    if (k == mpz_class(1)){return *this;}
    Point result(crv);
    Point base = *this;
    mpz_class exponent = k;

    while (exponent > 0) {
        if (exponent % 2 == 1) {
            result = result + base;
        }
        base = base + base;
        exponent /= 2;
    }

    return result;
}

void Point::set_order(mpz_class ord) {
    this->order = ord;
}

mpz_class Point::get_order() const{
    return this->order;
}



/*mpz_class Point::calculate_order() const {
    if (this->isInfinity()) {
        return mpz_class(1);  
    }

    if (this->get_y() == 0) {
        return mpz_class(2);  
    }
    mpz_class n(static_cast<long>(this->crv->points.size()));
    mpz_class m = sqrt(n) + 1;
    Point temp = *this * m;
    Point qpoint(temp.get_x(), -temp.get_y(), temp.get_curve());
    std::vector<Point> p_points, q_points;
    p_points.push_back(*this);
    q_points.push_back(qpoint);
    for (mpz_class i = 1; i < m; ++i) {
    size_t idx = i.get_ui(); 
    p_points.push_back(p_points[idx - 1] + *this);
    q_points.push_back(q_points[idx - 1] + qpoint);
    }

    for (mpz_class i = 0; i < m; ++i) {
        for (mpz_class j = 0; j < m; ++j) {
            if (p_points[i.get_ui()] == q_points[j.get_ui()]) {
                return i * m + j + 1;
            }
        }
    }
    return -1;
}*/

/*mpz_class Point::calculate_order() const {
    // Если точка на бесконечности, её порядок равен 1
    if (this->isInfinity())
        return mpz_class(1);
    
    // Если y == 0, то при удвоении получаем нейтральный элемент
    if (this->get_y() == 0)
        return mpz_class(2);
    
    mpz_class N = static_cast<unsigned long>(this->crv->points.size());
    
    mpz_class sqrt_N;
    mpz_sqrt(sqrt_N.get_mpz_t(), N.get_mpz_t());
    mpz_class m;
    if (sqrt_N * sqrt_N < N)
        m = sqrt_N + 1;
    else
        m = sqrt_N;
    
    auto point_to_string = [this](const Point& pt) -> std::string {
        if (pt.isInfinity())
            return "inf";
        return pt.get_x().get_str() + "," + pt.get_y().get_str();
    };
    
    // Шаг 1. Вычисляем "маленькие шаги": для j = 0, 1, ..., m-1 вычисляем jP
    // и сохраняем в map: ключ – строковое представление точки, значение – j.
    std::map<std::string, mpz_class> babySteps;
    for (mpz_class j = 0; j < m; j++) {
        Point baby = (*this) * j;  // 0*P, 1*P, 2*P, ..., (m-1)*P
        std::string key = point_to_string(baby);
        // Сохраним только первый (наименьший) индекс j для данного представления точки
        if (babySteps.find(key) == babySteps.end()) {
            babySteps[key] = j;
        }
    }
    
    Point factor = (*this) * m;
    for (mpz_class i = 1; i < m; i++) {
        Point giant = factor * i;
        Point neg_giant = giant.isInfinity()
                            ? giant
                            : Point(giant.get_x(), mod(-giant.get_y(), this->crv->get_p()), this->crv);
        
        std::string key = point_to_string(neg_giant);
        if (babySteps.find(key) != babySteps.end()) {
            mpz_class j = babySteps[key];
            mpz_class order = i * m + j;
            if (order > 0)
                return order;
        }
    }
    
    // Если ничего не найдено – выполняем запасной перебор (хотя такой случай встречается редко)
    mpz_class order = 1;
    Point current = *this;
    while (!current.isInfinity()) {
        order++;
        current = current + *this;
        if (order > N)
            break;
    }
    return order;
}*/

std::string Point::to_string() const {
    if (isInfinity()) return "inf";
    return x.get_str() + "," + y.get_str();
}

void Point::calculate_order() {
    if (this->isInfinity()) {
        this->set_order(mpz_class(1));
        return; 
    }
    if (this->get_y() == 0) {
        this->set_order(mpz_class(2));
        return;
    }

    mpz_class N = crv->get_group_order();

    // Вычисляем m = ceil(sqrt(N))
    mpz_class m;
    mpz_sqrt(m.get_mpz_t(), N.get_mpz_t());
    if (m * m < N) {
        m = m + 1;
    }

    std::vector<Point> babySteps;
    for (int i = 0; i < m; ++i) {
        babySteps.push_back(*this * i);
    }

    Point M = *this * m;

    std::vector<Point> giantSteps;
    for (int j = 0; j < m; ++j) {
        // Вычисляем j * M
        Point temp = M * j;
        if (!temp.isInfinity()) {
            mpz_class negY = crv->get_p() - temp.get_y();
            temp = Point(temp.get_x(), negY, temp.get_curve(), false);
        }
        giantSteps.push_back(temp);
    }

    bool found = false;
    mpz_class orderCandidate;
    for (int j = 0; j < m && !found; ++j) {
        for (int i = 0; i < m; ++i) {
            if (i == 0 && j == 0) continue; // пропускаем тривиальное совпадение Infinity == Infinity
            if (babySteps[i] == giantSteps[j]) {
                orderCandidate = i + m * j;
                found = true;
                break;
            }
        }
    }

    if (found) {
        this->set_order(orderCandidate);
    } else {
        for (mpz_class i = N; i > mpz_class(0); --i){
            if ((*this * i).isInfinity()) {
                this->set_order(N);
            }
        }
    }
}

