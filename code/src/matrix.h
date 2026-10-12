#ifndef MATRIX_H
#define MATRIX_H
#include <stddef.h>
double *matrix_alloc(size_t rows, size_t stride);
void matrix_free(double *p);
void matrix_mult_ptr_reg(size_t n, const double *a, size_t sa,
                         const double *b, size_t sb, double *c, size_t sc);
int matrix_mult_base_stride(size_t n, const double *a, size_t sa,
                            const double *b, size_t sb, double *c, size_t sc);
int strassen_stride(size_t n, const double *a, size_t sa, const double *b,
                    size_t sb, double *c, size_t sc, size_t cutoff);
int matrix_has_avx(void);
int matrix_compare(size_t n, const double *a, size_t sa, const double *b,
                    size_t sb, double atol, double rtol, double *max_error);
#endif
