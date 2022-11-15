#include <cufft.h>  
#include <cuComplex.h>


typedef float fp;
typedef cuComplex fpComplex;
typedef bfpNumFloat bfpNum;

typedef struct {
    unsigned short sign;   // 1 bit
    unsigned int exp;   // 8 bit
    int mant;   // 23 bit
} bfpNumFloat;

typedef struct{
    float original;
    unsigned short sign;   // 1 bit
    short exp;   // 11 bit
    long mant;   // 52 bit
} bfpNumDouble;

typedef struct {
    unsigned int common_exp;
    vector<int> M;
} bfpBlock;
