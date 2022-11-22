#pragma once
#include <cufft.h>  
#include <cuComplex.h>
#include <vector>


#define FLOAT_LEN_BITSIZE 32
#define FLOAT_EXP_BITSIZE 8
#define FLOAT_MANT_BITSIZE 23
#define DOUBLE_EXP_BITSIZE 11
#define DOUBLE_MANT_BITSIZE 52
#define MAX_BITSIZE 100

using namespace std;

typedef float fp;
typedef cuComplex fpComplex;

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
    vector<unsigned int> sign;
    vector<int> M;
} bfpBlock;
