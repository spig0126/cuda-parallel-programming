/* 
STRUCTURES
- BFP in float & double
- size of exponent and mantissa

FUNCTIONS
- 
- conversion to floating-point to bfp
    => divide floating-point into mantassa, exponent, sign
- put N FP numbers in one block
    => by class/list/pointer/structure
- choose common exponent
- mantissa alignment
- fixed-point operations on blocks

BFP structure - mantassa, exponent, sign,and block
*/

#include "data.cuh"
#include <string.h>

#define FLOAT_LEN_BITSIZE 32
#define FLOAT_EXP_BITSIZE 8
#define FLOAT_MANT_BITSIZE 23
#define DOUBLE_EXP_BITSIZE 11
#define DOUBLE_MANT_BITSIZE 52
#define MAX_BITSIZE 100

#define N 5

#define SIGN_EXTRACT_NUM -2147483648 
// #define MANT_EXTRACT_NUM 8388607
// #define EXP_EXTRACT_NUM 2139095040

/* 
short 2 byte
unsigned short 2byte
int 4 byte
unsign int 4 byte
float 4 byte
double 8 byte
*/
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

typedef bfpNumFloat bfpNum;
typedef float fp;

typedef struct {
    unsigned int common_exp;
    bfpNumFloat M[N];
} bfpBlock;


void printBit_bfpNumFloat(bfpNumFloat b, bool nextLine);
void printBit_uint(unsigned int num, int len, bool reverse);
float bfpNumFloat_to_float(bfpNumFloat b);
bfpNumFloat float_to_bfpNumFloat(float f);

/* bit representation operations */
void printBit_bfpNumFloat(bfpNumFloat b, bool nextLine){
    // printf("---------Bit representation of %f--------------\n", b.original);
    printBit_uint(b.sign, 1, false);
    printf(" ");
    printBit_uint(b.exp, FLOAT_EXP_BITSIZE, false);
    printf(" ");
    printBit_uint(b.mant, FLOAT_MANT_BITSIZE, false);

    if(nextLine){
        printf("\n");

    }
}

void printBit_uint(unsigned int num, int len, bool reverse){
    char out[MAX_BITSIZE] = "";

    for(int i=0; i<len; i++, num>>=1){
        if(num & 1){
            strcat(out, "1");
        }
        else{
            strcat(out, "0");
        }
    }

    if(reverse){
        for(int i=0; i<len; i++){
            printf("%c", out[i]);
        }
        
    }
    else{
        for(int i=len-1; i>=0; i--){
            printf("%c", out[i]);
        }
    }
}

void printBit_float(float f){
    printBit_bfpNumFloat(float_to_bfpNumFloat(f), true);
}

/* prints out signed integer with 4bits grouped together*/
void printBit_sint(int num, bool newLine){   
    int temp = num;
    char out[MAX_BITSIZE] = "";
    for(int i=0; i<32; i++, temp>>=1){
        if(i%4 == 0){
            strcat(out, " ");
        }
        if(temp&1){
            strcat(out, "1");
        }
        else{
            strcat(out, "0");
        }
    }
    
    for(int i=39; i>=0; i--){
         printf("%c", out[i]);
    }
    printf(" (0x%.8x)", num);

    if(newLine){
        printf("\n");
    }
}

void print_bfpNumFloat(bfpNumFloat b){
    printf("---------BFP component----------\n");
    printf("sign: %d\n", b.sign);
    printf("exp: %d\n", b.exp);
    printf("mant: %d\n", b.mant);
}

void print_bfpBlock(bfpBlock block){
    printf("----------BFP block ------------------\n");
    printf("N: %d\ncommon_exp: (%d)  ", N, block.common_exp);
    printBit_uint(block.common_exp, FLOAT_EXP_BITSIZE, false);
    printf("\n");
    for(int i=0; i<N; i++){
        printf("%d: ", i);
        printBit_bfpNumFloat(block.M[i], true);
    }
}




/* type conversions */
float bfpNumFloat_to_float(bfpNumFloat b){
    unsigned int t = 0;
    int sign = b.sign & 1;
    int exp = b.exp & 0x000000ff;
    int mant = b.mant & 0x007fffff;
    t ^= sign;
    t <<= FLOAT_EXP_BITSIZE;
    t ^= exp;
    t <<= FLOAT_MANT_BITSIZE;
    t ^= b.sign? ~mant + 1: mant;

    float* res = reinterpret_cast<float*>(&t);
    return *res;
}

bfpNumFloat float_to_bfpNumFloat(float f){
    bfpNumFloat res = {0};

    unsigned int* t = reinterpret_cast<unsigned int*>(&f);  //cast to integer for bitwise operations

    int mant_extract_num = 8388607;
    int exp_extract_num = 255;
    res.mant = *t & mant_extract_num;
    *t >>= FLOAT_MANT_BITSIZE;

    res.exp = *t & exp_extract_num;
    *t >>= FLOAT_EXP_BITSIZE;
    
    res.sign = *t & 1;
    // if(res.sign){
    //     res.mant = ~res.mant + 1;
    // }

    return res;
}


/* block formatting */
bfpBlock createBfpBlock(fp* X){
    bfpBlock block;

    //transform floating-point numbers into bfp format
    for(int i=0; i<N; i++){
        block.M[i] = float_to_bfpNumFloat(X[i]);
    } 

    //find common exponent
    unsigned int max_exp = 0;
    for(int i=0; i<N; i++){
        if(block.M[i].exp > max_exp){
            max_exp = block.M[i].exp;
        }
    }
    block.common_exp = max_exp;


    //align mantissas
    for(int i=0; i<N; i++){
        // printf("\nbefore: ");
        // printBit_uint(block.M[i].mant, 32, false);

        block.M[i].mant = (block.M[i].mant ^ 0x00800000) >> (max_exp - block.M[i].exp);

        // printf("\nafter:  ");
        // printBit_uint(block.M[i].mant, 32, false);
    }

    return block;
}

/* arithmetic operations */
bfpNumFloat add_f(bfpBlock block, bfpNumFloat a, bfpNumFloat b){
    bfpNumFloat res = {0};
    res.sign = 0;
    res.exp = block.common_exp;

    //2. change to 2's complement
    printf("\n\n1. after conversion to 2's complement\n");
    printf("\ta: ");
    printBit_sint(a.mant, false);
    printf(" => ");
    if(a.sign){
        a.mant = ~a.mant + 1;
    }
    printBit_sint(a.mant, true);
    printf("\tb: ");
    printBit_sint(b.mant, false);
    printf(" => ");
    if(b.sign){
        b.mant = ~b.mant + 1;
    }
    printBit_sint(b.mant, true);
    printf("\n");
    

    //2. add mantissas
    res.mant = a.mant + b.mant;
    printf("2. add mantissas\n\t");
    printBit_sint(res.mant, true);
    printf("\n");

    //3. convert to signed magnitude if negative

    if((res.mant & 0x80000000) == 0x80000000){    //if MSB is 1(sign = 1)
        res.mant = ~res.mant + 1;
        res.sign = 1;
    }
    printf("3. convert to signed magnitude if negative\n");
    printf("\t");
    printBit_sint(res.mant, true);
    printf("\n");

    //4. normalization
    int temp = res.mant >> FLOAT_MANT_BITSIZE + 1;
    if(temp & 1){
        res.mant >>= 1;
        res.exp += 1;
        temp >= 1;
    }
    printf("4. normalization\n");
    printf("\t");
    printBit_sint(res.mant, true);
    printf("\n");

    // printf("before mant: \n");
    // printf("\ta: ");
    // printBit_uint(block.M[a_i], 32, false);
    // printf("\n");
    // printf("\tb: ");
    // printBit_uint(block.M[b_i], 32, false);
    // printf("\n");

    // printf("after mant: \n");
    // printBit_uint(res.mant, 32, false);
    // printf("\n");


    //5. remove implicit 1
    printf("5. remove implicit 1\n\t");
    res.mant -= 0x00800000;
    printBit_uint(res.mant, 32, false);
    printf("\n");
    return res;
}

bfpNumFloat sub_f(bfpNumFloat a, bfpNumFloat b){
    bfpNum res = {0};
    if(a.exp != b.exp){
        printf("ADD ERROR: the two bfpNum don't have the same exponent\n");
        print_bfpNumFloat(a);
        print_bfpNumFloat(b);
        exit(0);
    }
    res.mant = a.mant - b.mant;
    res.sign = res.mant^SIGN_EXTRACT_NUM? 1: 0;
    res.exp = a.exp;
    return res;
}

bfpNumFloat mult_f(bfpNumFloat a, bfpNumFloat b){
}