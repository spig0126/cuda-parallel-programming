/* 
arranged modifications
- change bfpBlock components so that only mantissas are stored, not bfpNumFloat
- 
*/

#include "data.cuh"
#include <string.h>
#include <float.h>


#define FLOAT_LEN_BITSIZE 32
#define FLOAT_EXP_BITSIZE 8
#define FLOAT_MANT_BITSIZE 23
#define DOUBLE_EXP_BITSIZE 11
#define DOUBLE_MANT_BITSIZE 52
#define MAX_BITSIZE 100

#define N 2

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

/* prints out unsigned integer with 4bits grouped together*/
void printBit_ulong(long long num, bool newLine){
    long long temp = num;
    char out[MAX_BITSIZE] = "";
    for(int i=0; i<64; i++, temp>>=1){
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
    
    for(int i=79; i>=0; i--){
         printf("%c", out[i]);
    }
    // printf(" (0x%.8x)", num);

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
    t ^= mant;

    float* res = reinterpret_cast<float*>(&t);
    return *res;
}

bfpNumFloat float_to_bfpNumFloat(float f){
    bfpNumFloat res = {0};

    //cast to integer for bitwise operations
    unsigned int* t = reinterpret_cast<unsigned int*>(&f);  

    //extract mantissa bits(23 bits) using bitwise AND operation
    int mant_extract_num = 0x007fffff;
    res.mant = (*t & mant_extract_num) ^ 0x00800000;    //add implicit 1 in mantissa
    *t >>= FLOAT_MANT_BITSIZE;

    //extract exponent bits(8 bits) using bitwise AND operation
    int exp_extract_num = 0x000000ff;
    res.exp = *t & exp_extract_num;
    *t >>= FLOAT_EXP_BITSIZE;
    
    //extract sign bit using bitwise AND operation
    res.sign = *t & 1;

    return res;
}


/* block formatting */
bfpBlock createBfpBlock(fp* X){
    bfpBlock block;

    //transform floating-point numbers into bfp format, ㅁ
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
        block.M[i].mant = block.M[i].mant >> (max_exp - block.M[i].exp);
    }

    return block;
}

/* arithmetic operations for 2 numbers */
bfpNumFloat add_f(bfpBlock block, bfpNumFloat a, bfpNumFloat b){
    bfpNumFloat res = {0};
    res.sign = 0;
    res.exp = block.common_exp;

    //1. Conversion to 2’s complement for negative mantissas
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
    if(res.mant & 0x01000000){  //when carry is 1
        res.mant >>= 1;
        res.exp += 1;
    }
    while((res.mant & 0x00800000) != 0x00800000){   //11.01 x 2^1 = 1.101 x 2^2
        res.mant <<= 1;
        res.exp -= 1;
    }
    printf("4. normalization\n");
    printf("\t");
    printBit_sint(res.mant, true);
    printf("\n");

    //5. remove implicit 1
    res.mant &= 0x007fffff;
    printf("5. remove implicit 1\n\t");
    printBit_sint(res.mant, true);

    return res;
}

bfpNumFloat sub_f(bfpBlock block, bfpNumFloat a, bfpNumFloat b){
    bfpNumFloat b_neg = b;
    b_neg.sign = 1;
    return add_f(block, a, b_neg);
}

bfpNumFloat mult_f(bfpBlock block, bfpNumFloat a, bfpNumFloat b){
    bfpNumFloat res = {(unsigned short)(a.sign^b.sign), block.common_exp + block.common_exp - 127, 0};

    //1. multiply mantissas
    // unsigned long temp =   => 23bit * 23 bit = 46 bit, 앞 23bit만 가져오기
    //saturation function: 32bit 중에서 가장 큰 숫자 표현하기(overflow 발생했을 때 사용) /C++ intrinsic function -> only used in integer operations(not FP)
    //approximation compution
    //overflow/underflow 그냥 수용하기 (negative )
    unsigned long long res_temp = (unsigned long long)a.mant * (unsigned long long)b.mant;
    printf("1. multiply mantissas\n\t");
    printBit_ulong(res_temp, true);
    printf("\n");

    //2. normalization
    while((res_temp & 0x0000800000000000) != 0x0000800000000000 && (res_temp != 0)){
        res_temp <<= 1;
        res.exp-=1;
    }
    printf("\n2. normalization\n");
    printf("\t");
    printBit_ulong(res_temp, true);
    printf("\n");

    //3 rounding
    unsigned short last_bit =  (unsigned short)((res_temp & 0x0000000001000000) >> 24);
    unsigned short ground_bit = (unsigned short)((res_temp & 0x0000000000800000) >> 23);
    unsigned short round_bit =  (unsigned short)((res_temp & 0x0000000000400000) >> 22);
    unsigned short sticky_bits = (unsigned short)(res_temp & 0x00000000003fffff);

    res_temp &= 0x0000ffffff000000; //truncate
    if(ground_bit == 1){
        if(round_bit == 0 && sticky_bits == 0){ //round to even
            if(last_bit == 1){
                res_temp += 0x0000000001000000;
            }
        }
        else{
            res_temp += 0x0000000001000000; //round up
        }
    }
    printf("\n3. rounding\n");
    printf("\t");
    printBit_ulong(res_temp, true);
    printf("\t");
    for(int i=0; i<10; i++){
        printf("-----");
    }
    printf("^");
    for(int i=0; i<5; i++){
        printf("------");
    }
    printf("\n\n");

    //4. normalization(carry)
    res.mant = res_temp >> 23;
    int carry = res.mant >> 24;
    while(carry > 0){
        res.mant >>= 1;
        carry >>= 1;
        res.exp += 1;
    }
    printf("4. normalization(carry)\n");
    printf("\t");
    printBit_sint(res.mant, true);
    printf("\n");

    //5. remove implicit 1
    res.mant &= 0x007fffff;
    printf("5. remove implicit 1\n\t");
    printBit_sint(res.mant, true);

    return res;
}

bfpNumFloat div_f(bfpBlock block, bfpNumFloat a, bfpNumFloat b){
    /* 
    - need to consider when divisor=0
    */
    bfpNumFloat res = {(unsigned short)(a.sign^b.sign), 127, 0};

    //1. divide mantissas
    unsigned long long a_temp = (unsigned long long)a.mant<<40;
    unsigned long long b_temp = (unsigned long long)b.mant;
    unsigned long long res_temp = a_temp / b_temp;

    printf("\n1. divide mantissas\n\t");
    printBit_ulong(res_temp, true);
    printf("\n");

    // 2. normalization
    while((res_temp & 0x0000010000000000) != 0x0000010000000000 && (res_temp != 0)){
        res_temp <<= 1;
        res.exp-=1;
    }
    printf("\n2. normalization\n");
    printf("\t");
    printBit_ulong(res_temp, true);
    printf("\n");

    //3. rouding to nearest even
    unsigned short last_bit =  (unsigned short)((res_temp & 0x0000000000020000) >> 17);
    unsigned short ground_bit = (unsigned short)((res_temp & 0x0000000000010000) >> 16);
    unsigned short round_bit =  (unsigned short)((res_temp & 0x0000000000008000) >> 15);
    unsigned short sticky_bits = (unsigned short)(res_temp & 0x0000000000007fff);

    printf("\n3. rounding\n");
    printf("\t");
    printBit_ulong(res_temp, true);
    printf("\t");
    for(int i=0; i<11; i++){
        printf("-----");
    }
    printf("--^-");
    for(int i=0; i<4; i++){
        printf("------");
    }
    res_temp &= 0x000001fffffe0000; //truncate
    if(ground_bit == 1){
        if(round_bit == 0 && sticky_bits == 0){ //round to even
            if(last_bit == 1){
                res_temp += 0x0000000000020000;
            }
        }
        else{
            res_temp += 0x0000000000020000; //round up
        }
    }

    printf("\n\t=>\n\t");
    printBit_ulong(res_temp, true);
    printf("\t");
    for(int i=0; i<11; i++){
        printf("-----");
    }
    printf("--^-");
    for(int i=0; i<4; i++){
        printf("------");
    }
    printf("\n\n");

    //4.normalization (Carry)
    res.mant = (int)(res_temp >> 17);
    int carry = res.mant >> 24;
    while(carry > 0){
        res.mant >>= 1;
        carry >>= 1;
        res.exp += 1;
    }
    printf("\n4. normalization(carry)\n");
    printf("\t");
    printBit_sint(res.mant, true);
    printf("\n");

    //5. remove implicit 1
    res.mant &= 0x007fffff;
    printf("5. remove implicit 1\n\t");
    printBit_sint(res.mant, true);

    return res;
}
/* arithmetic operations for entire block */
bfpNumFloat add_block_f(bfpBlock block){
    bfpNumFloat M[N];
    bfpNumFloat res = {0, block.common_exp, 0};

    // copy block's mantissa list to M[N]
    for(int i=0; i<N; i++){
        M[i] = block.M[i];
    }
    // 1. conversion to 2's complement for negative mantissas
    for(int i=0; i<N; i++){
        if(M[i].sign){
            M[i].mant = ~M[i].mant + 1;
        }
    }

    //2. add mantissas
    for(int i=0; i<N; i++){
        res.mant += M[i].mant;
    }

    //3. convert to sign magnitude if addition result is negative
    if((res.mant & 0x80000000) == 0x80000000){
        res.mant = ~res.mant + 1;
        res.sign = 1;
    }

    //4. normalization
    int carry = res.mant >> 6;
    while(carry > 0){
        res.mant >>= 1;
        carry >>= 1;
        res.exp += 1;
    }
    while((res.mant & 0x00800000) != 0x00800000){
        res.mant <<= 1;
        res.exp -= 1;
    }

    //5. remove implicit 1
    res.mant &= 0xff7fffff;

    return res;
}


