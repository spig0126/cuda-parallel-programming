/* basic complex number operations */
#include <math.h>
#include <cuComplex.h>

__device__ __host__ cuDoubleComplex cal_euler(double x){
    cuDoubleComplex res;
    res.x = cos(x);
    res.y = sin(x);
    return res;
}
