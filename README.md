```markdown
# Реализация ГОСТ 34.10-2018 и Стрибог

Этот проект предоставляет реализацию алгоритмов цифровой подписи ГОСТ 34.10-2018 и хеш-функции Стрибог (GOST R 34.11-2012) на языке C. Проект включает операции с эллиптическими кривыми, хеширование данных и работу с большими числами с использованием библиотеки GMP.

## Содержание
- [Требования](#требования)
- [Сборка](#сборка)
- [Структура проекта](#структура-проекта)
- [Использование](#использование)
- [Примеры](#примеры)
- [Лицензия](#лицензия)

## Требования
- Компилятор C (например, `gcc`).
- Библиотека [GMP](https://gmplib.org/) для работы с большими числами.
- Стандартные библиотеки C.

## Сборка
1. Установите GMP, если она отсутствует:
   ```bash
   sudo apt-get install libgmp-dev  # Для Debian/Ubuntu
   ```
2. Скомпилируйте проект:
   ```bash
   gcc main.c src/hash/*.c src/sign/*.c src/ec/*.c -o gost3410.exe -lgmp
   ```
   Замените `main.c` на ваш файл с тестовым кодом.

## Структура проекта
- **ec_point.c/h**: Операции с точками эллиптических кривых (сложение, умножение, инициализация).
- **stribog.c/h**: Реализация хеш-функции Стрибог (поддержка 256 и 512-битных хешей).
- **stribog_data.h**: Константы и S-блоки для Стрибога.
- **gost3410.c/h**: Функции для создания и проверки цифровых подписей по ГОСТ 34.10-2018.
- **types.h**: Определения типов данных (u8, u16 и т.д.).

## Использование
### Хеширование данных
Пример использования хеш-функции Стрибог:
```c
#include "stribog.h"

struct stribog_ctx_t ctx;
u8 data[] = {0x01, 0x02, 0x03};
init(&ctx, HASH512); // Выбор размера хеша (HASH256 или HASH512)
stribog(&ctx, data, sizeof(data));
// Результат в ctx.h
```

### Подпись и проверка сообщения
```c
#include "gost3410.h"

mpz_t d, q, p, a, r, s;
EC_Point P, Q;
// Инициализация параметров кривой, ключей и точки P
gost3410_sign(r, s, message, message_len, d, q, p, a, &P);
int isValid = gost3410_verify(message, message_len, r, s, &Q, q, p, a, &P);
```

## Примеры
Пример работы с подписью:
```c
// Параметры эллиптической кривой (примерные значения)
mpz_init_set_str(p, "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F", 16);
mpz_init_set_str(a, "0", 16);
mpz_init_set_str(q, "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141", 16);

// Закрытый ключ d и базовая точка P
mpz_init_set_str(d, "1E240", 16);
ec_point_init(&P);
mpz_set_str(P.x, "79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798", 16);
mpz_set_str(P.y, "483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8", 16);

// Подпись сообщения
unsigned char msg[] = "Hello, GOST!";
gost3410_sign(r, s, msg, sizeof(msg)-1, d, q, p, a, &P);

// Проверка подписи
EC_Point Q;
ec_point_mul(&Q, d, &P, p, a); // Q = d * P
int isValid = gost3410_verify(msg, sizeof(msg)-1, r, s, &Q, q, p, a, &P);
printf("Подпись %s\n", isValid ? "верна" : "неверна");
```
