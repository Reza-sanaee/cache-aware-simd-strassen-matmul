#define _POSIX_C_SOURCE 200112L
#include "matrix.h"
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#ifdef _WIN32
#include <malloc.h>
#endif
#if defined(__GNUC__) && (defined(__x86_64__) || defined(__i386__)) && !defined(MATRIX_DISABLE_AVX)
#include <immintrin.h>
#define HAVE_AVX_KERNEL 1
#endif

double *matrix_alloc(size_t rows, size_t stride) {
    if (!rows || !stride || rows > SIZE_MAX / stride / sizeof(double)) return NULL;
    size_t bytes = rows * stride * sizeof(double);
#ifdef _WIN32
    return (double *)_aligned_malloc(bytes, 64);
#else
    void *p = NULL;
    return posix_memalign(&p, 64, bytes) == 0 ? p : NULL;
#endif
}
void matrix_free(double *p) {
#ifdef _WIN32
    _aligned_free(p);
#else
    free(p);
#endif
}
int matrix_has_avx(void) {
#ifdef HAVE_AVX_KERNEL
    __builtin_cpu_init();
    return __builtin_cpu_supports("avx") != 0;
#else
    return 0;
#endif
}
void matrix_mult_ptr_reg(size_t n, const double *a, size_t sa,
                         const double *b, size_t sb, double *c, size_t sc) {
    for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < n; ++j) {
            double sum = 0.0;
            for (size_t k = 0; k < n; ++k) sum += a[i*sa+k] * b[k*sb+j];
            c[i*sc+j] = sum;
        }
}
#ifdef HAVE_AVX_KERNEL
__attribute__((target("avx")))
static void leaf_avx(size_t n, const double *a, size_t sa,
                      const double *bt, double *c, size_t sc) {
    /* Sixteen vector accumulators evaluate a 4 x 4 output tile. */
    for (size_t i = 0; i < n; i += 4)
        for (size_t j = 0; j < n; j += 4) {
            size_t nr = n-i < 4 ? n-i : 4, nc = n-j < 4 ? n-j : 4;
            __m256d acc[4][4];
            for (size_t r=0;r<nr;++r)
                for(size_t q=0;q<nc;++q) acc[r][q]=_mm256_setzero_pd();
            size_t k=0;
            for (; k+4<=n; k+=4) {
                __m256d av[4], bv[4];
                for(size_t r=0;r<nr;++r) av[r]=_mm256_loadu_pd(a+(i+r)*sa+k);
                for(size_t q=0;q<nc;++q) bv[q]=_mm256_loadu_pd(bt+(j+q)*n+k);
                for(size_t r=0;r<nr;++r)
                    for(size_t q=0;q<nc;++q)
                        acc[r][q]=_mm256_add_pd(acc[r][q],_mm256_mul_pd(av[r],bv[q]));
            }
            for(size_t r=0;r<nr;++r)
                for(size_t q=0;q<nc;++q) {
                    double lanes[4]; _mm256_storeu_pd(lanes,acc[r][q]);
                    double sum=lanes[0]+lanes[1]+lanes[2]+lanes[3];
                    for(size_t tail=k;tail<n;++tail) sum+=a[(i+r)*sa+tail]*bt[(j+q)*n+tail];
                    c[(i+r)*sc+j+q]=sum;
                }
        }
}
#endif
int matrix_mult_base_stride(size_t n, const double *a, size_t sa,
                            const double *b, size_t sb, double *c, size_t sc) {
    if(!n || sa<n || sb<n || sc<n) return 0;
    double *bt=matrix_alloc(n,n);if(!bt) return 0;
    for(size_t i=0;i<n;++i)
        for(size_t j=0;j<n;++j) bt[j*n+i]=b[i*sb+j];
#ifdef HAVE_AVX_KERNEL
    if(matrix_has_avx()) leaf_avx(n,a,sa,bt,c,sc);
    else
#endif
    for(size_t i=0;i<n;++i)
        for(size_t j=0;j<n;++j) {
            double sum=0.0;
            for(size_t k=0;k<n;++k) sum+=a[i*sa+k]*bt[j*n+k];
            c[i*sc+j]=sum;
        }
    matrix_free(bt);return 1;
}
static void combine(size_t n,const double *a,size_t sa,const double *b,size_t sb,
                      double *c,int subtract) {
    for(size_t i=0;i<n;++i)
        for(size_t j=0;j<n;++j) c[i*n+j]=a[i*sa+j]+(subtract ? -b[i*sb+j] : b[i*sb+j]);
}
int strassen_stride(size_t n,const double *a,size_t sa,const double *b,size_t sb,
                    double *c,size_t sc,size_t cutoff) {
    if(!n || sa<n || sb<n || sc<n || !cutoff) return 0;
    if(n<=cutoff || n%2) return matrix_mult_base_stride(n,a,sa,b,sb,c,sc);
    size_t m=n/2;double *p[9]={0};int ok=0;
    for(size_t i=0;i<9;++i) if(!(p[i]=matrix_alloc(m,m))) goto cleanup;
    const double *a11=a,*a12=a+m,*a21=a+m*sa,*a22=a+m*sa+m;
    const double *b11=b,*b12=b+m,*b21=b+m*sb,*b22=b+m*sb+m;
    double *t1=p[7],*t2=p[8];
#define PRODUCT(index,x,sx,y,sy) do { if(!strassen_stride(m,x,sx,y,sy,p[index],m,cutoff)) goto cleanup; } while(0)
    combine(m,a11,sa,a22,sa,t1,0);combine(m,b11,sb,b22,sb,t2,0);PRODUCT(0,t1,m,t2,m);
    combine(m,a21,sa,a22,sa,t1,0);PRODUCT(1,t1,m,b11,sb);
    combine(m,b12,sb,b22,sb,t1,1);PRODUCT(2,a11,sa,t1,m);
    combine(m,b21,sb,b11,sb,t1,1);PRODUCT(3,a22,sa,t1,m);
    combine(m,a11,sa,a12,sa,t1,0);PRODUCT(4,t1,m,b22,sb);
    combine(m,a21,sa,a11,sa,t1,1);combine(m,b11,sb,b12,sb,t2,0);PRODUCT(5,t1,m,t2,m);
    combine(m,a12,sa,a22,sa,t1,1);combine(m,b21,sb,b22,sb,t2,0);PRODUCT(6,t1,m,t2,m);
#undef PRODUCT
    for(size_t i=0;i<m;++i)
        for(size_t j=0;j<m;++j) {
            size_t k=i*m+j,idx=i*sc+j;
            c[idx]=p[0][k]+p[3][k]-p[4][k]+p[6][k];
            c[idx+m]=p[2][k]+p[4][k];
            c[idx+m*sc]=p[1][k]+p[3][k];
            c[idx+m*sc+m]=p[0][k]-p[1][k]+p[2][k]+p[5][k];
        }
    ok=1;
cleanup:
    for(size_t i=0;i<9;++i) matrix_free(p[i]);
    return ok;
}
int matrix_compare(size_t n,const double *a,size_t sa,const double *b,size_t sb,
                    double atol,double rtol,double *max_error) {
    int ok=1;*max_error=0.0;
    for(size_t i=0;i<n;++i)
        for(size_t j=0;j<n;++j) {
            double x=a[i*sa+j],y=b[i*sb+j],error=fabs(x-y);
            if(!isfinite(x)||!isfinite(y)) {ok=0;*max_error=INFINITY;}
            else {
                if(error>*max_error)*max_error=error;
                if(error>atol+rtol*fmax(fabs(x),fabs(y)))ok=0;
            }
        }
    return ok;
}
