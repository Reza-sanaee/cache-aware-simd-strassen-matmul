#define _POSIX_C_SOURCE 200809L
#include "matrix.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <errno.h>
#ifdef _WIN32
#include <windows.h>
#else
#include <time.h>
#endif
static double seconds(void) {
#ifdef _WIN32
    LARGE_INTEGER t,f;QueryPerformanceCounter(&t);QueryPerformanceFrequency(&f);
    return (double)t.QuadPart/(double)f.QuadPart;
#else
    struct timespec t;clock_gettime(CLOCK_MONOTONIC,&t);return t.tv_sec+t.tv_nsec*1e-9;
#endif
}
static uint32_t next_random(uint32_t *s) {*s^=*s<<13;*s^=*s>>17;*s^=*s<<5;return *s;}
static int parse(const char *s,size_t *value) {
    if(!*s || *s=='-')return 0;
    char *end;errno=0;unsigned long long n=strtoull(s,&end,10);
    if(errno || *end || n<1 || n>8192)return 0;
    *value=(size_t)n;return 1;
}
int main(int argc,char **argv) {
    size_t n=128,cutoff=512,seed=42;
    if(argc>4 || (argc>1&&!parse(argv[1],&n)) || (argc>2&&!parse(argv[2],&cutoff)) || (argc>3&&!parse(argv[3],&seed))) {
        fprintf(stderr,"Usage: matrix_benchmark [dimension 1..8192] [cutoff 1..8192] [seed 1..8192]\n");return 2;
    }
    double *a=matrix_alloc(n,n),*b=matrix_alloc(n,n),*ref=matrix_alloc(n,n),*opt=matrix_alloc(n,n);
    if(!a||!b||!ref||!opt){fprintf(stderr,"Memory allocation failed\n");matrix_free(a);matrix_free(b);matrix_free(ref);matrix_free(opt);return 2;}
    uint32_t state=(uint32_t)seed;
    for(size_t i=0;i<n*n;++i) {
        a[i]=(double)(1+next_random(&state)%12345)/(1+next_random(&state)%6789);
        b[i]=(double)(1+next_random(&state)%12345)/(1+next_random(&state)%6789);
    }
    double t=seconds();matrix_mult_ptr_reg(n,a,n,b,n,ref,n);double reference=seconds()-t;
    t=seconds();int allocated=strassen_stride(n,a,n,b,n,opt,n,cutoff);double optimized=seconds()-t;
    double error=0.0;int pass=allocated&&matrix_compare(n,ref,n,opt,n,1e-8,1e-10,&error);
    puts("dimension,cutoff,seed,avx,reference_seconds,optimized_seconds,speedup,max_absolute_error,verified");
    printf("%zu,%zu,%zu,%d,%.9f,%.9f,%.6f,%.12g,%d\n",n,cutoff,seed,matrix_has_avx(),reference,optimized,optimized>0?reference/optimized:0,error,pass);
    matrix_free(a);matrix_free(b);matrix_free(ref);matrix_free(opt);return pass?0:1;
}
