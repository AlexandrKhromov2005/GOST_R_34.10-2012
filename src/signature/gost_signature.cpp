#include "gost_signature.hpp"
#include "../ec/algorithms_for_primes.hpp"
#include <fstream>
#include <random>

GOSTSignature::GOSTSignature(Curve* crv, const Point& basePoint, const mpz_class& subgroupOrder)
    : curve(crv),
      P(basePoint),  // Явно используем конструктор копирования
      q(subgroupOrder) 
{
    if (!curve->is_on_curve(P)) {
        throw std::invalid_argument("Точка не на кривой");
    }
}

// Генерация ключевой пары
void GOSTSignature::generateKeyPair() {
    std::random_device rd;
    gmp_randstate_t state;
    gmp_randinit_default(state);
    gmp_randseed_ui(state, rd());

    mpz_class max = q - 1;
    mpz_urandomm(d.get_mpz_t(), state, max.get_mpz_t());
    d += 1; // 0 < d < q

    Q = P * d; // Q = dP
}

// Формирование подписи
std::pair<mpz_class, mpz_class> GOSTSignature::sign(const u8* message, size_t messageLen) {
    // Хеширование сообщения
    u8 hash[OUTPUT_SIZE_512];
    struct stribog_ctx_t ctx;
    init(&ctx, HASH512);
    stribog(&ctx, const_cast<u8*>(message), messageLen);
    memcpy(hash, ctx.h, OUTPUT_SIZE_512);

    // Преобразование хеша в число e
    mpz_class e = bytesToMPZ(hash, OUTPUT_SIZE_512) % q;
    if (e == 0) e = 1;

    // Генерация случайного k
    mpz_class k;
    do {
        k = generateRandomNumber(q);
    } while (k == 0);

    // Вычисление точки C = kP
    Point C = P * k;
    mpz_class r = C.get_x() % q;
    if (r == 0) {
        throw std::runtime_error("r == 0, повторите подпись");
    }

    // Вычисление s = (r*d + k*e) mod q
    mpz_class s = (r * d + k * e) % q;
    if (s == 0) {
        throw std::runtime_error("s == 0, повторите подпись");
    }

    return {r, s};
}

// Проверка подписи
bool GOSTSignature::verify(const u8* message, size_t messageLen, const std::pair<mpz_class, mpz_class>& signature) {
    mpz_class r = signature.first;
    mpz_class s = signature.second;

    // Проверка 0 < r, s < q
    if (r <= 0 || r >= q || s <= 0 || s >= q) return false;

    // Хеширование сообщения
    u8 hash[OUTPUT_SIZE_512];
    struct stribog_ctx_t ctx;
    init(&ctx, HASH512);
    stribog(&ctx, const_cast<u8*>(message), messageLen);
    memcpy(hash, ctx.h, OUTPUT_SIZE_512);

    // Преобразование хеша в число e
    mpz_class e = bytesToMPZ(hash, OUTPUT_SIZE_512) % q;
    if (e == 0) e = 1;

    // Вычисление v = e^(-1) mod q
    mpz_class v = mod_inverse(e, q);

    // Вычисление z1 = s*v mod q, z2 = -r*v mod q
    mpz_class z1 = (s * v) % q;
    mpz_class z2 = (-r * v) % q;

    // Вычисление точки C' = z1*P + z2*Q
    Point C1 = P * z1;
    Point C2 = Q * z2;
    Point C_prime = C1 + C2;

    // Проверка r == x_C' mod q
    return (r == C_prime.get_x() % q);
}

mpz_class GOSTSignature::bytesToMPZ(const u8* bytes, size_t len) {
    mpz_class num;
    mpz_import(num.get_mpz_t(), len, 1, 1, 0, 0, bytes);
    return num;
}

mpz_class GOSTSignature::generateRandomNumber(const mpz_class& max) {
    gmp_randstate_t state;
    gmp_randinit_default(state);
    gmp_randseed_ui(state, std::random_device{}());
    
    mpz_class num;
    mpz_urandomm(num.get_mpz_t(), state, max.get_mpz_t());
    return num;
}