/* basic complex number operations */
#include <math.h>
#include "data.cuh"

__device__ __host__ fpComplex cal_euler(fp x){
    fpComplex res;
    res.x = cos(x);
    res.y = sin(x);
    return res;
}
