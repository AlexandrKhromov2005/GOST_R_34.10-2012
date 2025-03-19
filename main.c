#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>
#include "src/sign/gost3410.h"
#include "src/ec/ec_point.h"
#include "src/hash/types.h"
#include "src/hash/stribog.h"

// Тестовые параметры из задания
#define D_STR "7A929ADE789BB9BE10ED359DD39A72C11B60961F49397EEE1D19CE9891EC3B28"
#define QX_STR "7F2B49E270DB6D90D8595BEC458B50C58585BA1D4E9B788F6689DBD8E56FD80B"
#define QY_STR "26F1B489D6701DD185C8413A977B3CBBAF64D1C593D26627DFFB101A87FF77DA"

// Параметры кривой из ГОСТ Р 34.10-2012 (пример)
#define P_STR "8000000000000000000000000000000000000000000000000000000000000431"
#define A_STR "7"
#define Q_STR "8000000000000000000000000000000150FE8A1892976154C59CFC193ACCF5B3"
#define PX_STR "2"
#define PY_STR "8E2A8A0E65147D4BD6316030E16D19C85C97F0A9CA267122B96ABBCEA7E8FC8"

// Функция для чтения содержимого файла в буфер
unsigned char* read_file(const char *filename, size_t *out_size) {
    FILE *fp = fopen(filename, "rb");
    if (!fp) {
        perror("Error opening file");
        exit(EXIT_FAILURE);
    }
    fseek(fp, 0, SEEK_END);
    long size = ftell(fp);
    rewind(fp);
    unsigned char *buffer = (unsigned char*)malloc(size);
    if (!buffer) {
        fclose(fp);
        fprintf(stderr, "Memory allocation error\n");
        exit(EXIT_FAILURE);
    }
    if (fread(buffer, 1, size, fp) != (size_t)size) {
        fclose(fp);
        free(buffer);
        fprintf(stderr, "Error opening file\n");
        exit(EXIT_FAILURE);
    }
    fclose(fp);
    *out_size = size;
    return buffer;
}

// Функция записи подписи в файл sig.txt (r и s в 16-ричном виде, по одной строке)
void write_signature(const mpz_t r, const mpz_t s, const char *filename) {
    FILE *fp = fopen(filename, "w");
    if (!fp) {
        perror("Error opening file to write signature");
        exit(EXIT_FAILURE);
    }
    // Записываем r и s в 16-ричном виде
    gmp_fprintf(fp, "%Zx\n%Zx\n", r, s);
    fclose(fp);
}

// Функция чтения подписи из файла sig.txt
void read_signature(mpz_t r, mpz_t s, const char *filename) {
    FILE *fp = fopen(filename, "r");
    if (!fp) {
        perror("Error opening file with signature");
        exit(EXIT_FAILURE);
    }
    char buf[1024];
    if (!fgets(buf, sizeof(buf), fp)) {
        fclose(fp);
        fprintf(stderr, "Error reading signature (r)\n");
        exit(EXIT_FAILURE);
    }
    // Убираем символ перевода строки
    buf[strcspn(buf, "\r\n")] = 0;
    if (mpz_set_str(r, buf, 16) != 0) {
        fclose(fp);
        fprintf(stderr, "Error parsing signature (r)\n");
        exit(EXIT_FAILURE);
    }
    if (!fgets(buf, sizeof(buf), fp)) {
        fclose(fp);
        fprintf(stderr, "Error reading signature (s)\n");
        exit(EXIT_FAILURE);
    }
    buf[strcspn(buf, "\r\n")] = 0;
    if (mpz_set_str(s, buf, 16) != 0) {
        fclose(fp);
        fprintf(stderr, "Error parsing signature (r)\n");
        exit(EXIT_FAILURE);
    }
    fclose(fp);
}

int main() {
    // Инициализация параметров кривой
    mpz_t p, a, q;
    mpz_inits(p, a, q, NULL);
    mpz_set_str(p, P_STR, 16);
    mpz_set_str(a, A_STR, 16);
    mpz_set_str(q, Q_STR, 16);

    // Базовая точка P
    EC_Point P;
    ec_point_init(&P);
    mpz_set_str(P.x, PX_STR, 16);
    mpz_set_str(P.y, PY_STR, 16);
    P.infinity = 0;

    // Закрытый ключ d
    mpz_t d;
    mpz_init_set_str(d, D_STR, 16);

    // Публичный ключ Q (заданный тестовыми параметрами)
    EC_Point Q;
    ec_point_init(&Q);
    mpz_set_str(Q.x, QX_STR, 16);
    mpz_set_str(Q.y, QY_STR, 16);
    Q.infinity = 0;

    int choice;
    printf("Choose operation\n");
    printf("1. Sign file.txt and save the signature in sig.txt\n");
    printf("2. Check signature (file: file.txt, signature: sig.txt)\n");
    printf("Your choice: ");
    if (scanf("%d", &choice) != 1) {
        fprintf(stderr, "Input error\n");
        exit(EXIT_FAILURE);
    }
    
    // Очищаем буфер ввода
    int c;
    while ((c = getchar()) != '\n' && c != EOF);

    if (choice == 1) {
        // Подпись файла
        size_t msg_len;
        unsigned char *msg = read_file("file.txt", &msg_len);
        printf("=== DEBUG: Reading file.txt, size = %zu bytes ===\n", msg_len);
        
        mpz_t r, s;
        mpz_inits(r, s, NULL);

        printf("=== DEBUG: Call gost3410_sign ===\n");
        gost3410_sign(r, s, msg, msg_len, d, q, p, a, &P);
        gmp_printf("Signature:\nr = %Zx\ns = %Zx\n", r, s);
        
        write_signature(r, s, "sig.txt");
        printf("Signature saved in sig.txt\n");
        
        free(msg);
        mpz_clears(r, s, NULL);
    } else if (choice == 2) {
        // Проверка подписи
        size_t msg_len;
        unsigned char *msg = read_file("file.txt", &msg_len);
        printf("=== DEBUG: Reading file.txt, size = %zu bytes ===\n", msg_len);
        
        mpz_t r, s;
        mpz_inits(r, s, NULL);
        read_signature(r, s, "sig.txt");
        gmp_printf("Read signature:\nr = %Zx\ns = %Zx\n", r, s);
        
        printf("=== DEBUG: Call gost3410_verify ===\n");
        int valid = gost3410_verify(msg, msg_len, r, s, &Q, q, p, a, &P);
        printf("Signature %s\n", valid ? "VALID" : "INVALID");
        
        free(msg);
        mpz_clears(r, s, NULL);
    } else {
        printf("Wrong choice\n");
    }

    // Очистка ресурсов
    mpz_clears(p, a, q, d, NULL);
    ec_point_clear(&P);
    ec_point_clear(&Q);
    
    return 0;
}
