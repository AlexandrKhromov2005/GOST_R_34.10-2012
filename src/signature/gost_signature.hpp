#ifndef GOST_SIGNATURE_HPP
#define GOST_SIGNATURE_HPP

#include "../ec/curve.hpp"
#include "../hash/stribog.h"
#include <gmpxx.h>

class GOSTSignature {
private:
    Curve* curve;        // Кривая с параметрами ГОСТ Р 34.10-2012
    Point P;             // Базовая точка
    mpz_class q;         // Порядок подгруппы
    mpz_class d;         // Закрытый ключ
    Point Q;             // Открытый ключ

    static mpz_class bytesToMPZ(const u8* bytes, size_t len);
    static mpz_class generateRandomNumber(const mpz_class& max);

public:
    // Инициализация параметров кривой и базовой точки
    GOSTSignature(Curve* crv, const Point& basePoint, const mpz_class& subgroupOrder);

    // Генерация ключевой пары (d, Q = dP)
    void generateKeyPair();

    // Подпись сообщения (возвращает пару (r, s))
    std::pair<mpz_class, mpz_class> sign(const u8* message, size_t messageLen);

    // Проверка подписи
    bool verify(const u8* message, size_t messageLen, const std::pair<mpz_class, mpz_class>& signature);

    // Сериализация/десериализация ключей
    void savePrivateKey(const char* filename);
    void loadPrivateKey(const char* filename);
    void savePublicKey(const char* filename);
    void loadPublicKey(const char* filename);
};

#endif