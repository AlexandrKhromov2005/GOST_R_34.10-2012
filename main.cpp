/*#include "hash/stribog.h"
#include "hash/stribog_data.h"
#include "hash/types.h"
#include "hash/stribog.c"*/


//test launch for display working stribog
/*int main(int argc, char *argv[]) {
	process_input(argc, argv);
	return EXIT_SUCCESS;
}*/

#include "src/signature/gost_signature.hpp"

int main() {
    // Инициализация параметров кривой (пример)
    mpz_class p("FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFD97", 16);
    mpz_class a("FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFD94", 16);
    mpz_class b("00000000000000000000000000000000000000000000000000000000000000a6", 16);
    mpz_class q("FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF6C611070995AD10045841B09B761B893", 16);
    
    Curve curve(a, b, p);
    Point P = curve.find_point_of_order(q); // Базовая точка

    // Создание объекта подписи
    GOSTSignature signer(&curve, P, q);

    // Генерация ключей
    signer.generateKeyPair();
    signer.savePrivateKey("private.key");
    signer.savePublicKey("public.key");

    // Подпись сообщения
    const char* message = "Hello, ГОСТ!";
    auto signature = signer.sign(reinterpret_cast<const u8*>(message), strlen(message));

    // Проверка подписи
    bool isValid = signer.verify(reinterpret_cast<const u8*>(message), strlen(message), signature);
    std::cout << "Подпись " << (isValid ? "валидна" : "невалидна") << std::endl;

    return 0;
}