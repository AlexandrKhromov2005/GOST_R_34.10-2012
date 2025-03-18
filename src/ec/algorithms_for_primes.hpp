#ifndef ALGORITHMS_FOR_PRIMES_HPP
#define ALGORITHMS_FOR_PRIMES_HPP

#include <gmpxx.h>
#include <gmp.h>
#include <vector>

mpz_class mod(const mpz_class& a, const mpz_class& p);
mpz_class pow_mod(const mpz_class& base, const mpz_class& exponent, const mpz_class& modulus);
mpz_class mod_sqrt(const mpz_class& a, const mpz_class& p);
mpz_class extended_gcd(const mpz_class& a, const mpz_class& b, mpz_class& x, mpz_class& y);
mpz_class mod_inverse(const mpz_class& a, const mpz_class& m);
int legendre_symbol(const mpz_class& a, const mpz_class& p);
std::vector<mpz_class> factorize(const mpz_class& n);

#endif // ALGORITHMS_FOR_PRIMES_HPP