/*#include <stdio.h>
#include <string.h>

#include "types.h"
#include "stribog.h"

void print_hash(struct stribog_ctx_t *ctx)
{
u8 i;

for (i = 0; i < (ctx->size ? OUTPUT_SIZE_512 : OUTPUT_SIZE_256); i++)
	printf("%02x", ctx->h[i]);

putchar('\n');
}

int main(int argc, char *argv[])
{
struct stribog_ctx_t ctx;

u8 message1[63] = {
	0x32, 0x31, 0x30, 0x39, 0x38, 0x37, 0x36, 0x35,
	0x34, 0x33, 0x32, 0x31, 0x30, 0x39, 0x38, 0x37,
	0x36, 0x35, 0x34, 0x33, 0x32, 0x31, 0x30, 0x39,
	0x38, 0x37, 0x36, 0x35, 0x34, 0x33, 0x32, 0x31,
	0x30, 0x39, 0x38, 0x37, 0x36, 0x35, 0x34, 0x33,
	0x32, 0x31, 0x30, 0x39, 0x38, 0x37, 0x36, 0x35,
	0x34, 0x33, 0x32, 0x31, 0x30, 0x39, 0x38, 0x37,
	0x36, 0x35, 0x34, 0x33, 0x32, 0x31, 0x30
};

u8 message2[72] = {
	0xfb, 0xe2, 0xe5, 0xf0, 0xee, 0xe3, 0xc8, 0x20,
	0xfb, 0xea, 0xfa, 0xeb, 0xef, 0x20, 0xff, 0xfb,
	0xf0, 0xe1, 0xe0, 0xf0, 0xf5, 0x20, 0xe0, 0xed,
	0x20, 0xe8, 0xec, 0xe0, 0xeb, 0xe5, 0xf0, 0xf2,
	0xf1, 0x20, 0xff, 0xf0, 0xee, 0xec, 0x20, 0xf1,
	0x20, 0xfa, 0xf2, 0xfe, 0xe5, 0xe2, 0x20, 0x2c,
	0xe8, 0xf6, 0xf3, 0xed, 0xe2, 0x20, 0xe8, 0xe6,
	0xee, 0xe1, 0xe8, 0xf0, 0xf2, 0xd1, 0x20, 0x2c,
	0xe8, 0xf0, 0xf2, 0xe5, 0xe2, 0x20, 0xe5, 0xd1
};

printf("Message 1:\n");

init(&ctx, HASH512);
stribog(&ctx, message1, sizeof(message1));

printf("H^{512}: ");
print_hash(&ctx);

init(&ctx, HASH256);
stribog(&ctx, message1, sizeof(message1));

printf("H^{256}: ");
print_hash(&ctx);

printf("\nMessage 2:\n");

init(&ctx, HASH512);
stribog(&ctx, message2, sizeof(message2));

printf("H^{512}: ");
print_hash(&ctx);

init(&ctx, HASH256);
stribog(&ctx, message2, sizeof(message2));

printf("H^{256}: ");
print_hash(&ctx);

return 0;
}*/


#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <errno.h>

#include "types.h"
#include "stribog.h"

#define HELP_MSG "Usage: %s [OPTIONS] [INPUT]\n"\
				"Options:\n"\
				"  -s <string>  Hash a string\n"\
				"  -f <file>    Hash a file\n"\
				"  -b <size>    Hash size (256 or 512, default: 512)\n"\
				"  -h           Show this help\n"

void write_hash_to_file(struct stribog_ctx_t *ctx, const char *filename) {
FILE *file = fopen(filename, "w");
if (!file) {
	fprintf(stderr, "Error opening file %s: %s\n", filename, strerror(errno));
	exit(EXIT_FAILURE);
}
u8 i;
u8 output_size = ctx->size ? OUTPUT_SIZE_512 : OUTPUT_SIZE_256;
for (i = 0; i < output_size; i++) {
	fprintf(file, "%02x", ctx->h[i]);
}
fprintf(file, "\n");
fclose(file);
}

void print_hash(struct stribog_ctx_t *ctx) {
	u8 i;
	u8 output_size = ctx->size ? OUTPUT_SIZE_512 : OUTPUT_SIZE_256;
	
	for (i = 0; i < output_size; i++)
		printf("%02x", ctx->h[i]);
	putchar('\n');
}

void hash_data(struct stribog_ctx_t *ctx, const u8 *data, size_t len) {
init(ctx, ctx->size ? HASH512 : HASH256);
stribog(ctx, (u8*)data, len);
write_hash_to_file(ctx, "hash.txt"); // Заменяем print_hash
}

void hash_file(struct stribog_ctx_t *ctx, const char *filename) {
	FILE *file = fopen(filename, "rb");
	if (!file) {
		fprintf(stderr, "Error opening file: %s\n", strerror(errno));
		exit(EXIT_FAILURE);
	}

	fseek(file, 0, SEEK_END);
	long file_size = ftell(file);
	rewind(file);

	u8 *buffer = (u8*)malloc(file_size);
	if (!buffer) {
		fclose(file);
		fprintf(stderr, "Memory allocation failed\n");
		exit(EXIT_FAILURE);
	}

	if (fread(buffer, 1, file_size, file) != file_size) {
		fclose(file);
		free(buffer);
		fprintf(stderr, "File read error\n");
		exit(EXIT_FAILURE);
	}

	init(ctx, ctx->size ? HASH512 : HASH256);
	stribog(ctx, buffer, file_size);
	write_hash_to_file(ctx, "hash.txt");

	free(buffer);
	fclose(file);
}

void process_input(int argc, char *argv[]) {
	struct stribog_ctx_t ctx;
	ctx.size = HASH512; // По умолчанию 512-битный хеш
	char *input_string = NULL;
	char *filename = NULL;

	// Разбор аргументов командной строки
	for (int i = 1; i < argc; i++) {
		if (strcmp(argv[i], "-s") == 0) {
			if (++i >= argc) {
				fprintf(stderr, "Missing string argument\n");
				exit(EXIT_FAILURE);
			}
			input_string = argv[i];
		} else if (strcmp(argv[i], "-f") == 0) {
			if (++i >= argc) {
				fprintf(stderr, "Missing filename argument\n");
				exit(EXIT_FAILURE);
			}
			filename = argv[i];
		} else if (strcmp(argv[i], "-b") == 0) {
			if (++i >= argc) {
				fprintf(stderr, "Missing hash size argument\n");
				exit(EXIT_FAILURE);
			}
			ctx.size = (atoi(argv[i]) == 256) ? HASH256 : HASH512;
		} else if (strcmp(argv[i], "-h") == 0) {
			printf(HELP_MSG, argv[0]);
			exit(EXIT_SUCCESS);
		} else {
			fprintf(stderr, "Unknown option: %s\n", argv[i]);
			exit(EXIT_FAILURE);
		}
	}

	if (input_string && filename) {
		fprintf(stderr, "Cannot process both string and file input\n");
		exit(EXIT_FAILURE);
	}

	if (input_string) {
		hash_data(&ctx, (u8*)input_string, strlen(input_string));
	} else if (filename) {
		hash_file(&ctx, filename);
	} else {
		u8 message1[63] = {
		0x32, 0x31, 0x30, 0x39, 0x38, 0x37, 0x36, 0x35,
		0x34, 0x33, 0x32, 0x31, 0x30, 0x39, 0x38, 0x37,
		0x36, 0x35, 0x34, 0x33, 0x32, 0x31, 0x30, 0x39,
		0x38, 0x37, 0x36, 0x35, 0x34, 0x33, 0x32, 0x31,
		0x30, 0x39, 0x38, 0x37, 0x36, 0x35, 0x34, 0x33,
		0x32, 0x31, 0x30, 0x39, 0x38, 0x37, 0x36, 0x35,
		0x34, 0x33, 0x32, 0x31, 0x30, 0x39, 0x38, 0x37,
		0x36, 0x35, 0x34, 0x33, 0x32, 0x31, 0x30
	};

	u8 message2[72] = {
		0xfb, 0xe2, 0xe5, 0xf0, 0xee, 0xe3, 0xc8, 0x20,
		0xfb, 0xea, 0xfa, 0xeb, 0xef, 0x20, 0xff, 0xfb,
		0xf0, 0xe1, 0xe0, 0xf0, 0xf5, 0x20, 0xe0, 0xed,
		0x20, 0xe8, 0xec, 0xe0, 0xeb, 0xe5, 0xf0, 0xf2,
		0xf1, 0x20, 0xff, 0xf0, 0xee, 0xec, 0x20, 0xf1,
		0x20, 0xfa, 0xf2, 0xfe, 0xe5, 0xe2, 0x20, 0x2c,
		0xe8, 0xf6, 0xf3, 0xed, 0xe2, 0x20, 0xe8, 0xe6,
		0xee, 0xe1, 0xe8, 0xf0, 0xf2, 0xd1, 0x20, 0x2c,
		0xe8, 0xf0, 0xf2, 0xe5, 0xe2, 0x20, 0xe5, 0xd1
	};

		printf("Original test vectors:\n");
		
		printf("Message 1:\n");
		init(&ctx, HASH512);
		stribog(&ctx, message1, sizeof(message1));
		printf("H^{512}: ");
		print_hash(&ctx);

		init(&ctx, HASH256);
		stribog(&ctx, message1, sizeof(message1));
		printf("H^{256}: ");
		print_hash(&ctx);

		printf("\nMessage 2:\n");
		init(&ctx, HASH512);
		stribog(&ctx, message2, sizeof(message2));
		printf("H^{512}: ");
		print_hash(&ctx);

		init(&ctx, HASH256);
		stribog(&ctx, message2, sizeof(message2));
		printf("H^{256}: ");
		print_hash(&ctx);
	}
}