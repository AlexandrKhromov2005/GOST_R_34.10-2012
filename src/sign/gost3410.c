#include "gost3410.h"
#include "../hash/stribog.h"    // Предполагается, что здесь объявлены init() и stribog()
#include "../hash/types.h"      // Определения u8, u64 и т.п.
#include "../ec/ec_point.h"
#include <gmp.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>

#define BLOCK_SIZE 64
#define GREEN   "\033[0;32m"
#define RED     "\033[0;31m"
#define RESET   "\033[0m"

/* Преобразование 64-байтного хэша в число mpz_t */
static void hash_to_mpz(mpz_t out, const unsigned char hash[BLOCK_SIZE]) {
    mpz_import(out, BLOCK_SIZE, 1, sizeof(unsigned char), 0, 0, hash);
}

/* Формирование подписи по ГОСТ 34.10–2018 */
void gost3410_sign(mpz_t r, mpz_t s,
                   const unsigned char *message, size_t message_len,
                   const mpz_t d, const mpz_t q,
                   const mpz_t p, const mpz_t a,
                   const EC_Point *P) {
    struct stribog_ctx_t ctx;
    init(&ctx, HASH512);
    stribog(&ctx, (u8 *)message, (u64)message_len);
    unsigned char hash[BLOCK_SIZE];
    memcpy(hash, ctx.h, BLOCK_SIZE);

    mpz_t a_value, e;
    mpz_inits(a_value, e, NULL);
    mpz_import(a_value, 64, 1, 1, 0, 0, hash);
    mpz_mod(e, a_value, q);
    if (mpz_cmp_ui(e, 0) == 0)
        mpz_set_ui(e, 1);

    //gmp_printf(GREEN "=== DEBUG: e = %Zx ===\n" RESET, e);

    gmp_randstate_t rand_state;
    gmp_randinit_default(rand_state);
    gmp_randseed_ui(rand_state, (unsigned long) time(NULL));

    EC_Point C;
    ec_point_init(&C);

    mpz_t k, rd, ke, temp;
    mpz_inits(k, rd, ke, temp, NULL);

    while (1) {
        mpz_urandomm(k, rand_state, q);
        if (mpz_cmp_ui(k, 0) == 0)
            continue;
        ec_point_mul(&C, k, P, p, a);
        if (C.infinity)
            continue;
        mpz_mod(r, C.x, q);
        if (mpz_cmp_ui(r, 0) == 0)
            continue;
        mpz_mul(rd, r, d);
        mpz_mul(ke, k, e);
        mpz_add(temp, rd, ke);
        mpz_mod(s, temp, q);
        if (mpz_cmp_ui(s, 0) == 0)
            continue;
        break;
    }

    //gmp_printf(GREEN "k = %Zx\n" RESET, k);
    //gmp_printf(GREEN "r = %Zx\n" RESET, r);
    //gmp_printf(GREEN "s = %Zx\n" RESET, s);

    ec_point_clear(&C);
    mpz_clears(a_value, e, k, rd, ke, temp, NULL);
    gmp_randclear(rand_state);
}

/* Проверка подписи по ГОСТ 34.10–2018 */
int gost3410_verify(const unsigned char *message, size_t message_len,
                    const mpz_t r, const mpz_t s,
                    const EC_Point *Q,
                    const mpz_t q, const mpz_t p, const mpz_t a,
                    const EC_Point *P) {
    if (mpz_cmp_ui(r, 0) <= 0 || mpz_cmp(r, q) >= 0 ||
        mpz_cmp_ui(s, 0) <= 0 || mpz_cmp(s, q) >= 0)
        return 0;

    struct stribog_ctx_t ctx;
    init(&ctx, HASH512);
    stribog(&ctx, (u8 *)message, (u64)message_len);
    unsigned char hash[BLOCK_SIZE];
    memcpy(hash, ctx.h, BLOCK_SIZE);

    //printf(RED "=== DEBUG: Hash of message ===\n" RESET);
    //for (int i = 0; i < BLOCK_SIZE; i++) printf(RED "%02x" RESET, hash[i]);
    //printf("\n\n");

    mpz_t a_value, e, v, z1, z2, temp;
    mpz_inits(a_value, e, v, z1, z2, temp, NULL);
    hash_to_mpz(a_value, hash);
    mpz_mod(e, a_value, q);
    if (mpz_cmp_ui(e, 0) == 0)
        mpz_set_ui(e, 1);

    //gmp_printf(RED "e = %Zx\n" RESET, e);

    if (mpz_invert(v, e, q) == 0) {
        //printf(RED "=== DEBUG: Inverse of e does not exist! ===\n" RESET);
        mpz_clears(a_value, e, v, z1, z2, temp, NULL);
        return 0;
    }

    mpz_mul(z1, s, v);
    mpz_mod(z1, z1, q);
    mpz_neg(z2, r);
    mpz_mul(z2, z2, v);
    mpz_mod(z2, z2, q);

    //gmp_printf(RED "z1 = %Zx\n" RESET, z1);
    //gmp_printf(RED "z2 = %Zx\n\n" RESET, z2);

    EC_Point R_point, temp1, temp2;
    ec_point_init(&R_point);
    ec_point_init(&temp1);
    ec_point_init(&temp2);
    ec_point_mul(&temp1, z1, P, p, a);
    ec_point_mul(&temp2, z2, Q, p, a);
    ec_point_add(&R_point, &temp1, &temp2, p, a);
    mpz_mod(temp, R_point.x, q);

    //gmp_printf(RED "R.x mod q = %Zx\n" RESET, temp);
    //gmp_printf(RED "r from signature = %Zx\n\n" RESET, r);

    int valid = (mpz_cmp(temp, r) == 0);
    //printf(RED "=== DEBUG: Signature VALID? %s ===\n" RESET, valid ? "YES" : "NO");

    ec_point_clear(&R_point);
    ec_point_clear(&temp1);
    ec_point_clear(&temp2);
    mpz_clears(a_value, e, v, z1, z2, temp, NULL);
    return valid;
}
