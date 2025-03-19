#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>
#include "src/sign/gost3410.h"
#include "src/ec/ec_point.h"
#include "src/hash/types.h"
#include "src/hash/stribog.h"

// Test parameters from the task
#define D_STR "7A929ADE789BB9BE10ED359DD39A72C11B60961F49397EEE1D19CE9891EC3B28"
// Known public key values for signature verification
#define QX_STR "7F2B49E270DB6D90D8595BEC458B50C58585BA1D4E9B788F6689DBD8E56FD80B"
#define QY_STR "26F1B489D6701DD185C8413A977B3CBBAF64D1C593D26627DFFB101A87FF77DA"

// Curve parameters from GOST R 34.10-2012 (example)
#define P_STR "8000000000000000000000000000000000000000000000000000000000000431"
#define A_STR "7"
#define Q_STR "8000000000000000000000000000000150FE8A1892976154C59CFC193ACCF5B3"
#define PX_STR "2"
#define PY_STR "8E2A8A0E65147D4BD6316030E16D19C85C97F0A9CA267122B96ABBCEA7E8FC8"

// Block size for Stribog
#define BLOCK_SIZE 64

/* Function to read the contents of a file into a buffer */
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
        fprintf(stderr, "Error reading file\n");
        exit(EXIT_FAILURE);
    }
    fclose(fp);
    *out_size = size;
    return buffer;
}

/* Function to write the signature to a file sig.txt (r and s in hexadecimal format, each on a separate line) */
void write_signature(const mpz_t r, const mpz_t s, const char *filename) {
    FILE *fp = fopen(filename, "w");
    if (!fp) {
        perror("Error opening file to write signature");
        exit(EXIT_FAILURE);
    }
    gmp_fprintf(fp, "%Zx\n%Zx\n", r, s);
    fclose(fp);
}

/* Function to read the signature from the file sig.txt */
void read_signature(mpz_t r, mpz_t s, const char *filename) {
    FILE *fp = fopen(filename, "r");
    if (!fp) {
        perror("Error opening signature file");
        exit(EXIT_FAILURE);
    }
    char buf[1024];
    if (!fgets(buf, sizeof(buf), fp)) {
        fclose(fp);
        fprintf(stderr, "Error reading signature (r)\n");
        exit(EXIT_FAILURE);
    }
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
        fprintf(stderr, "Error parsing signature (s)\n");
        exit(EXIT_FAILURE);
    }
    fclose(fp);
}

/* Main program */
int main() {
    // Initialize curve parameters
    mpz_t p, a, q;
    mpz_inits(p, a, q, NULL);
    mpz_set_str(p, P_STR, 16);
    mpz_set_str(a, A_STR, 16);
    mpz_set_str(q, Q_STR, 16);

    // Base point P
    EC_Point P;
    ec_point_init(&P);
    mpz_set_str(P.x, PX_STR, 16);
    mpz_set_str(P.y, PY_STR, 16);
    P.infinity = 0;

    // Private key d
    mpz_t d;
    mpz_init_set_str(d, D_STR, 16);

    // Calculate public key Q = d * P (coordinates are printed for debugging purposes)
    EC_Point Q_calc;
    ec_point_init(&Q_calc);
    ec_point_mul(&Q_calc, d, &P, p, a);
    gmp_printf("Calculated coordinates of Q:\nQ.x = %Zx\nQ.y = %Zx\n\n", Q_calc.x, Q_calc.y);

    // For verification, use the known public key value
    EC_Point Q;
    ec_point_init(&Q);
    mpz_set_str(Q.x, QX_STR, 16);
    mpz_set_str(Q.y, QY_STR, 16);
    Q.infinity = 0;

    int choice;
    printf("Choose an operation:\n");
    printf("1. Sign file file.txt and save the signature to sig.txt\n");
    printf("2. Verify signature (file: file.txt, signature: sig.txt)\n");
    printf("Your choice: ");
    if (scanf("%d", &choice) != 1) {
        fprintf(stderr, "Input error\n");
        exit(EXIT_FAILURE);
    }
    
    int c;
    while ((c = getchar()) != '\n' && c != EOF);

    if (choice == 1) {
        // Read the message from file
        size_t msg_len;
        unsigned char *msg = read_file("file.txt", &msg_len);
        printf("=== DEBUG: Read file.txt, size = %zu bytes ===\n", msg_len);
        
        mpz_t r, s;
        mpz_inits(r, s, NULL);

        printf("=== DEBUG: Calling gost3410_sign ===\n");
        gost3410_sign(r, s, msg, msg_len, d, q, p, a, &P);
        gmp_printf("Signature:\nr = %Zx\ns = %Zx\n", r, s);
        
        // Save the signature to file
        write_signature(r, s, "sig.txt");
        printf("Signature saved to sig.txt\n");
        
        free(msg);
        mpz_clears(r, s, NULL);
    } else if (choice == 2) {
        // Read the message from file
        size_t msg_len;
        unsigned char *msg = read_file("file.txt", &msg_len);
        printf("=== DEBUG: Read file.txt, size = %zu bytes ===\n", msg_len);
        
        mpz_t r, s;
        mpz_inits(r, s, NULL);
        read_signature(r, s, "sig.txt");
        gmp_printf("Read signature:\nr = %Zx\ns = %Zx\n", r, s);
        
        printf("=== DEBUG: Calling gost3410_verify ===\n");
        int valid = gost3410_verify(msg, msg_len, r, s, &Q, q, p, a, &P);
        printf("Signature is %s\n", valid ? "VALID" : "INVALID");
        
        free(msg);
        mpz_clears(r, s, NULL);
    } else {
        printf("Invalid choice\n");
    }

    // Clean up resources
    mpz_clears(p, a, q, d, NULL);
    ec_point_clear(&P);
    ec_point_clear(&Q);
    ec_point_clear(&Q_calc);
    
    return 0;
}
