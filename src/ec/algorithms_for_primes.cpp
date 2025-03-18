#include "algorithms_for_primes.hpp"
#include <algorithm>

mpz_class mod(const mpz_class &a, const mpz_class &p) {
    return (((a % p) + p) % p);
}

mpz_class mod_inverse(const mpz_class& a, const mpz_class& m) {
    mpz_class gcd, x, y;
    gcd = extended_gcd(a, m, x, y); // Добавлено присваивание gcd
    if (gcd != 1) throw std::runtime_error("Inverse does not exist");
    return mod(x, m);
}

mpz_class pow_mod(const mpz_class& base, const mpz_class& exponent, const mpz_class& modulus) {
    if (modulus == 1) return 0; // Ноль для modulus = 1
    mpz_class result = 1;
    mpz_class b = base % modulus;
    mpz_class e = exponent;

    while (e > 0) {
        if (e % 2 == 1) {
            result = (result * b) % modulus;
        }
        b = (b * b) % modulus;
        e /= 2;
    }
    return result;
}

mpz_class mod_sqrt(const mpz_class& a, const mpz_class& p) {
    if (a == 0) return 0; // Корень из нуля

    // Используем нашу функцию legendre_symbol вместо mpz_legendre
    if (legendre_symbol(a, p) != 1) return 0; // Не квадратичный вычет

    // Случай p = 2
    if (p == 2) return a;

    mpz_class q = p - 1;
    mpz_class s = 0;
    while (q % 2 == 0) {
        q /= 2;
        s++;
    }

    // Поиск невычета z с помощью нашей legendre_symbol
    mpz_class z = 2;
    while (legendre_symbol(z, p) != -1) {
        z++;
    }

    mpz_class c = pow_mod(z, q, p);
    mpz_class x = pow_mod(a, (q + 1) / 2, p);
    mpz_class t = pow_mod(a, q, p);
    mpz_class m = s;

    while (t != 1) {
        mpz_class tt = t;
        mpz_class i = 0;
        // Находим наименьшее i: t^(2^i) ≡ 1 mod p
        while (tt != 1 && i < m) {
            tt = pow_mod(tt, 2, p);
            i++;
        }

        // Вычисляем 2^(m-i-1) через GMP-функции
        mpz_class exponent = 1;
        mpz_class shift = m - i - 1;
        mpz_mul_2exp(exponent.get_mpz_t(), exponent.get_mpz_t(), mpz_get_ui(shift.get_mpz_t()));
        
        mpz_class b = pow_mod(c, exponent, p);
        x = (x * b) % p;
        t = (t * b % p * b) % p;
        c = pow_mod(b, 2, p);
        m = i;
    }

    return x;
}

mpz_class extended_gcd(const mpz_class& a, const mpz_class& b, mpz_class& x, mpz_class& y) {
    if (b == 0) {
        x = 1;
        y = 0;
        return a;
    }
    mpz_class q = a / b;
    mpz_class r = a % b;
    mpz_class x1, y1;
    mpz_class gcd = extended_gcd(b, r, x1, y1);
    x = y1;
    y = x1 - q * y1;
    return gcd;
}

bool invert(mpz_class& result, const mpz_class& a, const mpz_class& m) {
    if (m <= 1) {
        return false; 
    }
    mpz_class gcd, x, y;
    gcd = extended_gcd(a, m, x, y);
    if (gcd != 1) {
        return false; 
    }

    x %= m;
    if (x < 0) {
        x += m;
    }
    result = x;
    return true;
}


int legendre_symbol(const mpz_class& a, const mpz_class& p) {
    // Обработка p = 2
    if (p == 2) {
        mpz_class tmp = a % 2;
        return (tmp == 0) ? 0 : 1;
    }

    // Проверка на четность p (символ Лежандра требует нечетного простого p)
    if (p % 2 == 0) {
        return 0; // Некорректный ввод
    }

    mpz_class tmp = a % p;
    if (tmp == 0) return 0;

    int result = 1;
    mpz_class n = tmp;
    mpz_class d = p;

    while (n != 0) {
        // Извлечение множителей 2 из n
        unsigned long t = mpz_scan1(n.get_mpz_t(), 0); // Найти степень 2
        if (t > 0) {
            mpz_class divisor;
            mpz_ui_pow_ui(divisor.get_mpz_t(), 2, t); // divisor = 2^t
            n /= divisor; // Явное деление вместо n >>= t

            // Проверка условия для множителя 2
            mpz_class mod8 = d % 8;
            if ((t % 2 == 1) && (mod8 == 3 || mod8 == 5)) {
                result *= -1;
            }
        }

        // Квадратичный закон взаимности
        mpz_class mod4_n = n % 4;
        mpz_class mod4_d = d % 4;
        if (mod4_d == 3 && mod4_n == 3) {
            result *= -1;
        }

        // Обмен n и d с взятием модуля
        std::swap(n, d);
        n %= d;
    }

    return (d == 1) ? result : 0;
}

// Вспомогательная функция для алгоритма Полларда-Ро
mpz_class pollards_rho(const mpz_class& n) {
    if (n == 1) return 1;
    if (n % 2 == 0) return 2;

    mpz_class x = 2, y = 2, d = 1;
    mpz_class c = rand() % (n-1) + 1;

    // Явное объявление функтора
    class Lambda {
        mpz_class c, n;
    public:
        Lambda(mpz_class c, mpz_class n) : c(c), n(n) {}
        mpz_class operator()(const mpz_class& x) const {
            return (x * x + c) % n;
        }
    } f(c, n);

    while (d == 1) {
        x = f(x);
        y = f(f(y));
        mpz_sub(d.get_mpz_t(), x.get_mpz_t(), y.get_mpz_t());
        mpz_gcd(d.get_mpz_t(), d.get_mpz_t(), n.get_mpz_t());
    }

    return (d == n) ? mpz_class(0) : d;
}

// Рекурсивная факторизация с использованием Полларда-Ро
void factorize_recursive(const mpz_class& n, std::vector<mpz_class>& factors) {
    if (n == 1) return;
    if (mpz_probab_prime_p(n.get_mpz_t(), 15) > 0) {
        factors.push_back(n);
        return;
    }

    mpz_class d = pollards_rho(n);
    factorize_recursive(d, factors);
    factorize_recursive(n / d, factors);
}

// Основная функция факторизации

#include <algorithm>

std::vector<mpz_class> factorize(const mpz_class& n) {
    std::vector<mpz_class> factors;
    if (n == 0) return factors;

    mpz_class num = n;

    // Обработка отрицательных чисел
    if (num < 0) {
        factors.push_back(-1);
        num = -num;
    }

    // Удаляем множители 2
    while (num % 2 == 0) {
        factors.push_back(2);
        num /= 2;
    }

    // Проверяем нечётные делители до sqrt(num)
    mpz_class i(3);
    mpz_class max_i;
    mpz_sqrt(max_i.get_mpz_t(), num.get_mpz_t());
    
    while (i <= max_i) {
        while (num % i == 0) {
            factors.push_back(i);
            num /= i;
            mpz_sqrt(max_i.get_mpz_t(), num.get_mpz_t()); // Обновляем max_i
        }
        i += 2;
    }

    // Обработка оставшегося числа
    if (num > 1) {
        if (mpz_probab_prime_p(num.get_mpz_t(), 15) > 0) {
            factors.push_back(num);
        } else {
            // Используем Полларда-Ро и рекурсивно факторизуем
            mpz_class d = pollards_rho(num);
            if (d == 0 || d == 1) {
                factors.push_back(num); // Не удалось разложить
            } else {
                auto sub_factors = factorize(d);
                factors.insert(factors.end(), sub_factors.begin(), sub_factors.end());
                auto sub_remaining = factorize(num / d);
                factors.insert(factors.end(), sub_remaining.begin(), sub_remaining.end());
            }
        }
    }

    // Сортируем множители
    std::sort(factors.begin(), factors.end());
    return factors;
}

