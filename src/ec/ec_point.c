#include "ec_point.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

/* Инициализация точки: выделяются памяти для координат и устанавливается флаг бесконечности */
void ec_point_init(EC_Point *P) {
    mpz_init(P->x);
    mpz_init(P->y);
    P->infinity = 1; // по умолчанию точка на бесконечности
}

/* Очистка mpz_t */
void ec_point_clear(EC_Point *P) {
    mpz_clear(P->x);
    mpz_clear(P->y);
}

/* Копирование точки Q в R */
void ec_point_copy(EC_Point *R, const EC_Point *Q) {
    mpz_set(R->x, Q->x);
    mpz_set(R->y, Q->y);
    R->infinity = Q->infinity;
}

/* Функция модульного вычитания: вычисляет rop = a mod m, результат неотрицательный */
void mod_mpz(mpz_t rop, const mpz_t a, const mpz_t m) {
    mpz_mod(rop, a, m);
    if(mpz_sgn(rop) < 0) {
        mpz_add(rop, rop, m);
    }
}

/* Функция сложения точек на эллиптической кривой
   Параметры:
     - R – результат (выходная точка)
     - P, Q – входные точки
     - p – модуль конечного поля
     - a – коэффициент кривой
*/
void ec_point_add(EC_Point *R, const EC_Point *P, const EC_Point *Q, const mpz_t p, const mpz_t a) {
    // Если одна из точек на бесконечности, R = другая точка.
    if(P->infinity) {
        ec_point_copy(R, Q);
        return;
    }
    if(Q->infinity) {
        ec_point_copy(R, P);
        return;
    }
    
    // Если x1 == x2 и (y1 + y2) mod p == 0, то P + Q = O (точка на бесконечности)
    mpz_t temp;
    mpz_init(temp);
    mpz_add(temp, P->y, Q->y);
    mod_mpz(temp, temp, p);
    if(mpz_cmp(P->x, Q->x) == 0 && mpz_cmp_ui(temp, 0) == 0) {
        R->infinity = 1;
        mpz_clear(temp);
        return;
    }
    
    mpz_t lambda, numerator, denominator, inv;
    mpz_inits(lambda, numerator, denominator, inv, NULL);
    
    if(mpz_cmp(P->x, Q->x) == 0 && mpz_cmp(P->y, Q->y) == 0) {
        // Удвоение точки: lambda = (3*x1^2 + a) * (2*y1)^{-1} mod p
        mpz_mul(numerator, P->x, P->x);       // x1^2
        mpz_mul_ui(numerator, numerator, 3);    // 3*x1^2
        mpz_add(numerator, numerator, a);       // 3*x1^2 + a
        
        mpz_mul_ui(denominator, P->y, 2);         // 2*y1
    } else {
        // Сложение точек: lambda = (y2 - y1) * (x2 - x1)^{-1} mod p
        mpz_sub(numerator, Q->y, P->y);
        mpz_sub(denominator, Q->x, P->x);
    }
    
    mod_mpz(numerator, numerator, p);
    mod_mpz(denominator, denominator, p);
    
    // Вычисляем обратный элемент к denominator по модулю p
    if(mpz_invert(inv, denominator, p) == 0) {
        // Обратного элемента не существует, ошибка (обычно не происходит, если p – простое)
        fprintf(stderr, "Error: Inverse element not possible\n");
        exit(EXIT_FAILURE);
    }
    mpz_mul(lambda, numerator, inv);
    mod_mpz(lambda, lambda, p);
    
    // x3 = lambda^2 - x1 - x2 mod p
    mpz_t x3, y3;
    mpz_inits(x3, y3, NULL);
    mpz_mul(x3, lambda, lambda);
    mpz_sub(x3, x3, P->x);
    mpz_sub(x3, x3, Q->x);
    mod_mpz(x3, x3, p);
    
    // y3 = lambda*(x1 - x3) - y1 mod p
    mpz_sub(y3, P->x, x3);
    mpz_mul(y3, lambda, y3);
    mpz_sub(y3, y3, P->y);
    mod_mpz(y3, y3, p);
    
    // Записываем результат в R
    mpz_set(R->x, x3);
    mpz_set(R->y, y3);
    R->infinity = 0;
    
    mpz_clears(temp, lambda, numerator, denominator, inv, x3, y3, NULL);
}

/* Скалярное умножение: вычисляем R = k * P методом «двоичного разложения».
   Параметры:
     - R – результат (инициализируется внутри)
     - k – скаляр (mpz_t)
     - P – исходная точка
     - p – модуль конечного поля
     - a – коэффициент кривой (используется при сложении)
*/
void ec_point_mul(EC_Point *R, const mpz_t k, const EC_Point *P, const mpz_t p, const mpz_t a) {
    // Инициализируем R как точку на бесконечности
    ec_point_init(R);
    R->infinity = 1;
    
    EC_Point temp, base;
    ec_point_init(&base);
    ec_point_copy(&base, P);
    
    mpz_t exp;
    mpz_init_set(exp, k);
    
    // Алгоритм «double and add»
    while(mpz_cmp_ui(exp, 0) > 0) {
        if(mpz_odd_p(exp)) {
            EC_Point sum;
            ec_point_init(&sum);
            ec_point_add(&sum, R, &base, p, a);
            ec_point_copy(R, &sum);
            ec_point_clear(&sum);
        }
        EC_Point dbl;
        ec_point_init(&dbl);
        ec_point_add(&dbl, &base, &base, p, a);
        ec_point_copy(&base, &dbl);
        ec_point_clear(&dbl);
        mpz_fdiv_q_2exp(exp, exp, 1);  // делим exp на 2
    }
    
    ec_point_clear(&base);
    mpz_clear(exp);
}
