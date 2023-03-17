#ifndef UTILITY_BFP_H
#define UTILITY_BFP_H

#include "bfp.cuh"

using namespace bfp;

inline bfpNum degrees_to_radians(bfpNum degrees)
{
    return degrees * (b_pi / b_180);
}

inline float random_float()
{
    // Returns a random real in [0,1).
    return rand() / (RAND_MAX + 1.0);
}

inline float random_float(float min, float max)
{
    // Returns a random real in [min,max).
    return min + (max - min) * random_float();
}

inline bfpNum random_num()
{
    return float_to_bfpNum(random_float());
}

inline bfpNum random_num(bfpNum min, bfpNum max)
{
    float f_min = bfpNum_to_float(min);
    float f_max = bfpNum_to_float(max);

    return float_to_bfpNum(random_float(f_min, f_max));
}

inline float clamp(float x, float min, float max)
{
    if (x < min)
        return min;
    if (x > max)
        return max;
    return x;
}

inline bfpNum clamp(bfpNum x, bfpNum min, bfpNum max)
{
    if (x < min)
        return min;
    if (x > max)
        return max;
    return x;
}

// inline int random_int(int min, int max)
// {
//     // Returns a random integer in [min,max].
//     return static_cast<int>(random_num(min, max + 1));
// }
#endif