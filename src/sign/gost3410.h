#ifndef GOST3410_H
#define GOST3410_H

#include <stddef.h>
#include <gmp.h>
#include "../ec/ec_point.h"

/*
 * Функция формирования цифровой подписи ГОСТ 34.10–2018.
 *
 * Параметры:
 *   r, s         – mpz_t для выходных значений подписи.
 *   message      – сообщение для подписи (массив байт) и его длина.
 *   d            – закрытый ключ (0 < d < q).
 *   q            – порядок подгруппы эллиптической кривой.
 *   p, a         – параметры кривой: модуль конечного поля и коэффициент a.
 *   P            – базовая точка кривой.
 */
void gost3410_sign(mpz_t r, mpz_t s,
                   const unsigned char *message, size_t message_len,
                   const mpz_t d, const mpz_t q,
                   const mpz_t p, const mpz_t a,
                   const EC_Point *P);

/*
 * Функция проверки цифровой подписи ГОСТ 34.10–2018.
 *
 * Параметры:
 *   message      – сообщение, которое подписывали.
 *   r, s         – компоненты подписи.
 *   Q            – публичный ключ (Q = d*P).
 *   q            – порядок подгруппы эллиптической кривой.
 *   p, a         – параметры кривой.
 *   P            – базовая точка кривой.
 *
 * Возвращает 1, если подпись корректна, 0 – иначе.
 */
int gost3410_verify(const unsigned char *message, size_t message_len,
                    const mpz_t r, const mpz_t s,
                    const EC_Point *Q,
                    const mpz_t q, const mpz_t p, const mpz_t a,
                    const EC_Point *P);

#endif // GOST3410_H
