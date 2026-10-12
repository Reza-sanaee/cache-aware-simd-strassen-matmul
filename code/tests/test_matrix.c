#include "matrix.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
int main(void) {
    size_t sizes[]={1,2,3,4,5,7,8,12,16,31,32,33,64};size_t cases=0;
    for(size_t x=0;x<sizeof(sizes)/sizeof(sizes[0]);++x) {
        size_t n=sizes[x],sa=n+3,sb=n+5,sc=n+7;
        double *a=matrix_alloc(n,sa),*b=matrix_alloc(n,sb),*r=matrix_alloc(n,sc),*o=matrix_alloc(n,sc);
        assert(a&&b&&r&&o);
        for(int kind=0;kind<4;++kind) {
            for(size_t i=0;i<n;++i)
                for(size_t j=0;j<n;++j) {
                    a[i*sa+j]=kind==0?0:kind==1?(i==j):((int)((i*13+j*7)%31)-15)/7.0;
                    b[i*sb+j]=kind==2?(i==j):((int)((i*5+j*11)%29)-14)/9.0;
                }
            for(size_t i=0;i<n*sc;++i)o[i]=1234567.0;
            matrix_mult_ptr_reg(n,a,sa,b,sb,r,sc);
            assert(strassen_stride(n,a,sa,b,sb,o,sc,4));
            double error;assert(matrix_compare(n,r,sc,o,sc,1e-10,1e-10,&error));
            for(size_t i=0;i<n;++i)for(size_t j=n;j<sc;++j)assert(o[i*sc+j]==1234567.0);
            ++cases;
        }
        matrix_free(a);matrix_free(b);matrix_free(r);matrix_free(o);
    }
    double z=0.0,nan=NAN,error;assert(!matrix_compare(1,&z,1,&nan,1,1e-10,1e-10,&error));
    assert(matrix_compare(1,&z,1,&z,1,1e-10,1e-10,&error));
    printf("PASS: %zu matrix cases; padded-stride guards; zero and nonfinite comparison; AVX=%d\n",cases,matrix_has_avx());
    return 0;
}
