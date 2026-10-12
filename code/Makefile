CC = gcc
CFLAGS = -O3 -std=c11 -Wall -Wextra -Werror -Isrc
LDLIBS = -lm
.PHONY: all test clean
all: bin/matrix_benchmark
bin:
	mkdir -p bin
bin/matrix_benchmark: src/matrix.c src/main.c src/matrix.h | bin
	$(CC) $(CFLAGS) src/matrix.c src/main.c -o $@ $(LDLIBS)
bin/test_matrix: src/matrix.c tests/test_matrix.c src/matrix.h | bin
	$(CC) $(CFLAGS) src/matrix.c tests/test_matrix.c -o $@ $(LDLIBS)
bin/test_matrix_scalar: src/matrix.c tests/test_matrix.c src/matrix.h | bin
	$(CC) $(CFLAGS) -DMATRIX_DISABLE_AVX src/matrix.c tests/test_matrix.c -o $@ $(LDLIBS)
test: bin/test_matrix bin/test_matrix_scalar
	./bin/test_matrix
	./bin/test_matrix_scalar
clean:
	rm -f bin/matrix_benchmark bin/test_matrix bin/test_matrix_scalar
