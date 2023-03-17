#pragma once
#include <vector>


#define FLOAT_LEN_BITSIZE 32
#define FLOAT_EXP_BITSIZE 8
#define FLOAT_MANT_BITSIZE 23
#define DOUBLE_EXP_BITSIZE 11
#define DOUBLE_MANT_BITSIZE 52
#define MAX_BITSIZE 100

#define BFP_MANT_BITSIZE 15
#define BFP_SIGNIFICAND_BITSIZE 23
#define BFP_EXP_BITSIZE 8
#define BFP_BIAS (unsigned int)(std::pow(2, BFP_EXP_BITSIZE - 1) - 1)

using namespace std;


typedef float fp;

typedef struct {
    unsigned short sign;  
    unsigned int exp;  
    int mant;  
} bfpNum;

typedef struct {
    unsigned int common_exp;
    vector<unsigned int> sign;
    vector<int> M;
} bfpBlock;
