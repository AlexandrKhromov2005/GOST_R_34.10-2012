#ifndef EC_POINT_H
#define EC_POINT_H

#include <gmp.h>

/* Структура для представления точки эллиптической кривой.
   Поле infinity принимает значение 1 для точки на бесконечности, 0 – для обычной точки. */
typedef struct {
    mpz_t x;
    mpz_t y;
    int infinity;
} EC_Point;

/* Инициализация точки (выделяются mpz_t для координат) */
void ec_point_init(EC_Point *P);

/* Очистка ресурсов, выделенных в точке */
void ec_point_clear(EC_Point *P);

/* Копирование точки: R = Q */
void ec_point_copy(EC_Point *R, const EC_Point *Q);

/* Функция для вычисления a mod m, результат положительный */
void mod_mpz(mpz_t rop, const mpz_t a, const mpz_t m);

/* Сложение двух точек P и Q на эллиптической кривой по модулю p с коэффициентом a.
   Если одна из точек является точкой на бесконечности, возвращается другая.
   Если P == -Q, возвращается точка на бесконечности. */
void ec_point_add(EC_Point *R, const EC_Point *P, const EC_Point *Q, const mpz_t p, const mpz_t a);

/* Скалярное умножение: вычисление R = k * P с использованием метода «двоичного разложения».
   Все операции выполняются по модулю p. */
void ec_point_mul(EC_Point *R, const mpz_t k, const EC_Point *P, const mpz_t p, const mpz_t a);

#endif // EC_POINT_H
