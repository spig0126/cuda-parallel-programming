#ifndef UTILITY_BFP_H
#define UTILITY_BFP_H

#include "bfp.cuh"

using namespace bfp;

inline bfpNum degrees_to_radians(bfpNum degrees)
{
    return degrees * (b_pi / b_180);
}

inline bfpNum random_num()
{
    return float_to_bfpNum(rand() / (RAND_MAX + 1.0));
}

inline bfpNum random_num(bfpNum min, bfpNum max)
{
    return min + (max - min) * random_num();
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