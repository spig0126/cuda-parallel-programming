/* basic complex number operations */
#include <math.h>

typedef struct Comp{
    double r, i;
} Comp;

__device__ __host__ Comp cal_euler(double x){
    Comp res;
    res.r = cos(x);
    res.i = sin(x);
    return res;
}

__device__ __host__ Comp comp_mult(Comp a, Comp b){
    Comp res;
    res.r = a.r * b.r - a.i * b.i;
    res.i = a.r * b.i + a.i * b.r;
    return res;
}

__device__ __host__ Comp comp_add(Comp a, Comp b){
    Comp res;   
    res.r = a.r + b.r;
    res.i = a.i + b.i;
    return res;
}

__device__ __host__ Comp comp_sub(Comp a, Comp b){
    Comp res;   
    res.r = a.r - b.r;
    res.i = a.i - b.i;
    return res;
}
