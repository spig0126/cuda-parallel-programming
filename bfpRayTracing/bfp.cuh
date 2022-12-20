#include <float.h>
#include <math.h>

#include "bfpStruct.cuh"
#include "print.cuh"
#include "vec3.h"

using namespace std;

int BFP_MANT_EXTRACT_BITS = (int)pow(2, BFP_MANT_BITSIZE) - 1; 
int BFP_EXP_EXTRACT_BITS = (int)pow(2, BFP_EXP_BITSIZE) - 1;
int FLOAT_MANT_EXTRACT_BITS = 0x7fffff;
int FLOAT_EXP_EXTRACT_BITS = 0x7f800000;
int FLOAT_SIGN_EXTRACT_BITS = 0x80000000;
long long int LONG_ALL_ONES = 0xffffffffffffffff;
int INT_ALL_ONES = 0xffffffff;

/* type conversions */
float bfpNum_to_float(bfpNum b){
    unsigned int t = 0;
    int sign, exp, mant;

    sign = (int)b.sign;
    exp = (int)b.exp;
    mant = (int)b.mant;

    if(FLOAT_EXP_BITSIZE < BFP_EXP_BITSIZE){  
        exp >>= (BFP_EXP_BITSIZE - FLOAT_EXP_BITSIZE);
    }
    else{
        exp <<= (FLOAT_EXP_BITSIZE - BFP_EXP_BITSIZE);
    }
    if(FLOAT_MANT_BITSIZE < BFP_MANT_BITSIZE){
        mant >>= (BFP_MANT_BITSIZE - FLOAT_MANT_BITSIZE);
    }
    else{
        mant <<= (FLOAT_MANT_BITSIZE - BFP_MANT_BITSIZE);
    }
    
    mant ^= 0x800000;   //remove implicit 1

    t ^= sign;
    t <<= FLOAT_EXP_BITSIZE;
    t ^= exp;
    t <<= FLOAT_MANT_BITSIZE;
    t ^= mant;

    float* res = reinterpret_cast<float*>(&t);
    return *res;
}

bfpNum float_to_bfpNum(float f){
    bfpNum res = {0};

    //cast to integer for bitwise operations
    unsigned int* t = reinterpret_cast<unsigned int*>(&f);  

    //extract mantissa bits(23 bits) using bitwise AND operation
    res.mant = (*t & FLOAT_MANT_EXTRACT_BITS) ^ 0x00800000;    //add implicit 1 in mantissa

    //set bfp's mantissa from float's extracted mantissa bits
    if(FLOAT_MANT_BITSIZE >= BFP_MANT_BITSIZE){    //if bfp mantissa is smaller in bit size than float mantissa
        res.mant >>= (FLOAT_MANT_BITSIZE - BFP_MANT_BITSIZE);
    }

    //extract exponent bits(8 bits) using bitwise AND operation
    res.exp = (*t & FLOAT_EXP_EXTRACT_BITS) >> FLOAT_MANT_BITSIZE;

    //set bfp's exponent from float's extracted exponent bits
    if(FLOAT_EXP_BITSIZE >= BFP_EXP_BITSIZE){   //if bfp exponent is smaller in bit size than float exponent
        res.exp >>= (FLOAT_EXP_BITSIZE - BFP_EXP_BITSIZE);
    }
    
    //extract sign bit using bitwise AND operation
    res.sign = (*t & FLOAT_SIGN_EXTRACT_BITS) >> (FLOAT_MANT_BITSIZE + FLOAT_EXP_BITSIZE);

    return res;
}

vector<bfpNum> color_to_bfpNum(color c){
    vector<bfpNum> c_formatted;

    for(int i=0; i<3; i++){
        c_formatted.push_back(float_to_bfpNum(c[i]));
    }

    return c_formatted;
}

color bfpNums_to_color(vector<bfpNum> v){
    float r = bfpNum_to_float(v[0]);
    float g = bfpNum_to_float(v[1]);
    float b = bfpNum_to_float(v[2]);

    return color(r, g, b);
}


/* block formatting */
bfpBlock createBfpBlock(vector<float> X){
    bfpBlock block;

    //transform floating-point numbers into bfp format
    vector<bfpNum> temp;
    for(int i=0; i<X.size(); i++){
        temp.push_back(float_to_bfpNum(X[i])); //adds implicit 1 in the process!
    } 

    //find common exponent
    unsigned int max_exp = 0;
    for(int i=0; i<X.size(); i++){
        if(temp[i].exp > max_exp){
            max_exp = temp[i].exp;
        }
    }
    block.common_exp = max_exp;

    //save aligned mantissas to block
    for(int i=0; i<temp.size(); i++){
        block.M.push_back(temp[i].mant >> (max_exp - temp[i].exp));
        block.sign.push_back(temp[i].sign);
    }

    return block;
}

bfpBlock createColorBfpBlock(vector<color> colors){
    bfpBlock block;

    //transform floating-point numbers into bfp format
    vector<bfpNum> temp;
    for(int i=0; i<colors.size(); i++){
        vector<bfpNum> c_formatted = color_to_bfpNum(colors[i]);  //adds implicit 1 in the process!
        for(int j=0; j<3; j++){
            temp.push_back(c_formatted[j]);
        }
    } 

    //find common exponent
    unsigned int max_exp = 0;
    for(int i=0; i<temp.size(); i++){
        if(temp[i].exp > max_exp){
            max_exp = temp[i].exp;
        }
    }

    //save common exponent
    block.common_exp = max_exp;

    //save signs and aligned mantissas to block
    for(int i=0; i<temp.size(); i++){
        block.sign.push_back(temp[i].sign);
        block.M.push_back(temp[i].mant >> (max_exp - temp[i].exp));
    }

    return block;
}


/* arithmetic operations for 2 numbers */
bfpNum div_f(bfpBlock block, bfpNum a, bfpNum b){
    bfpNum res = {(unsigned short)(a.sign^b.sign), 127, 0};

    //1. divide mantissas
    unsigned long long a_temp = (unsigned long long)a.mant<<40;
    unsigned long long b_temp = (unsigned long long)b.mant;
    unsigned long long res_temp = a_temp / b_temp;

    // 2. normalization
    while((res_temp & 0x0000010000000000) != 0x0000010000000000 && (res_temp != 0)){
        res_temp <<= 1;
        res.exp-=1;
    }

    //3. rouding to nearest even
    unsigned short last_bit =  (unsigned short)((res_temp & 0x0000000000020000) >> 17);
    unsigned short ground_bit = (unsigned short)((res_temp & 0x0000000000010000) >> 16);
    unsigned short round_bit =  (unsigned short)((res_temp & 0x0000000000008000) >> 15);
    unsigned short sticky_bits = (unsigned short)(res_temp & 0x0000000000007fff);
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

    //4.normalization (Carry)
    res.mant = (int)(res_temp >> 17);
    int carry = res.mant >> 24;
    while(carry > 0){
        res.mant >>= 1;
        carry >>= 1;
        res.exp += 1;
    }

    //5. remove implicit 1
    res.mant &= 0x007fffff;

    return res;
}

/* arithmetic operations for entire block */
//float block
float add_bfpBlock(bfpBlock block){
    bfpNum res = {0, block.common_exp, 0};
    long long int res_mant_temp = 0;

    //1. converision to 2's complment for negative mantissas
    for(int i=0; i<block.sign.size(); i++){
        if(block.sign[i]){
            block.M[i] = ~block.M[i] + 1;
        }
    }

    //2. add mantissas
    for(int i=0; i<block.M.size(); i++){       
        res_mant_temp += (long long int)block.M[i];
        // printf("add ");
        // printBit_sint(block.M[i], false);
        // printf(" => ");
        // printBit_sint(res.mant, true);

        /* option 1 */
        // if((res.mant & 0xff000000) > 0){ //when carry is 1: 11.01 x 2^1 = 1.101 x 2^2
        //     res.mant >>= 1;
        //     res.exp += 1;
        //     printf("\tCARRY! : ");
        //     printBit_mant(res.mant, true);
        // }
    }

    //3. convert to  signed magnitude if negative
    if((res_mant_temp & 0x8000000000000000) == 0x8000000000000000){    //if MSB is 1(sign = 1)
        res_mant_temp = ~res_mant_temp + 1;
        res.sign = 1;
    }

    //4. normalization
    /* option 2 */
    while((res_mant_temp >> (BFP_MANT_BITSIZE + 1)) > 0){ //when carry is 1: 11.01 x 2^1 = 1.101 x 2^2
        res_mant_temp >>= 1;
        res.exp += 1;
        // printf("\tCARRY! : ");
        // printBit_bfpNum_exp(res.exp, false);
        // printBit_ulong(res_mant_temp, true);
    }
    res.mant = (int)res_mant_temp;
    int implicit_1 = (int)(pow(2, BFP_MANT_BITSIZE));
    while((res.mant & implicit_1) != implicit_1){   //0.101 x 2^2 = 1.01 x 2^1
        res.mant <<= 1;
        res.exp -= 1;
    }

    return bfpNum_to_float(res);
}

float mult_bfpBlock(bfpBlock block){
    bfpNum res = {0, (unsigned int)(block.common_exp * block.M.size() - ((int)pow(2, BFP_EXP_BITSIZE - 1) - 1) * (block.M.size()- 1)), 0};

    //1. assign sign
    for(int i=0; i<block.sign.size(); i++){
        res.sign ^= block.sign[i];
    }

    unsigned long long res_mant_temp = (unsigned long long)block.M[0];

    for(int i=1; i<block.M.size(); i++){
        //2. multiply mantissas
        res_mant_temp *= (unsigned long long) block.M[i];

        //3 rounding
        int temp = (int)pow(2, BFP_MANT_BITSIZE); 
        unsigned short last_bit =  (unsigned short)((res_mant_temp & temp) >> BFP_MANT_BITSIZE);
        temp >>= 1;
        unsigned short ground_bit = (unsigned short)((res_mant_temp & temp) >> BFP_MANT_BITSIZE - 1);
        temp >>= 1;
        unsigned short round_bit =  (unsigned short)((res_mant_temp & temp) >> BFP_MANT_BITSIZE - 2);
        temp -= 1;
        unsigned short sticky_bits = (unsigned short)(res_mant_temp & temp);

        int lsb = (int)pow(2, BFP_MANT_BITSIZE); 
        if(ground_bit == 1){
            if(round_bit == 0 && sticky_bits == 0){ //round to even
                if(last_bit == 1){
                    res_mant_temp += lsb;
                }
            }
            else{
                res_mant_temp += lsb; //round up
            }
        }

        res_mant_temp >>= BFP_MANT_BITSIZE;
    }

    //3. normalization
    while((res_mant_temp >> (BFP_MANT_BITSIZE + 1)) > 0){ //when carry is 1: 11.01 x 2^1 = 1.101 x 2^2
        res_mant_temp >>= 1;
        res.exp += 1;
        // printf("\tCARRY! : ");
        // printBit_bfpNum_exp(res.exp, false);
        // printBit_ulong(res_mant_temp, true);
    }    
    res.mant = (int)res_mant_temp;
    int implicit_1 = (int)(pow(2, BFP_MANT_BITSIZE));
    while((res.mant & implicit_1) != implicit_1){   //0.101 x 2^2 = 1.01 x 2^1
        res.mant <<= 1;
        res.exp -= 1;
    }


    //6. remove implicit 1
    res.mant = (int)res_mant_temp & 0x007fffff;

    
    return bfpNum_to_float(res);
}

//color block
color add_color_bfpBlock(bfpBlock block){
    vector<bfpNum> res(3, {0, block.common_exp, 0});
    vector<long long int> res_mant_temp(3, 0);

    //1. converision to 2's complment for negative mantissas
    for(int i=0; i<block.sign.size(); i++){
        if(block.sign[i]){
            block.M[i] = ~block.M[i] + 1;
        }
    }

    //2. add mantissas
    for(int i=0; i<block.M.size(); i+=3){
        res_mant_temp[0] += (long long int)block.M[i];
        res_mant_temp[1] += (long long int)block.M[i + 1];
        res_mant_temp[2] += (long long int)block.M[i + 2];

        // printBit_sint(res[0].mant, true);
        // printBit_sint(res[1].mant, true);
        // printBit_sint(res[2].mant, true);
    }

    //3. convert to  signed magnitude if negative
    for(int i=0; i<3; i++){
        if((res_mant_temp[i] & 0x8000000000000000) == 0x8000000000000000){    //if MSB is 1(sign = 1)
            res_mant_temp[i] = ~res_mant_temp[i] + 1;
            res[i].sign = 1;
        }
    }

    //4. normalization
    for(int i=0; i<3; i++){
        while((res_mant_temp[i] >> (BFP_MANT_BITSIZE + 1)) > 0){ //when carry is 1
            res_mant_temp[i] >>= 1;
            res[i].exp += 1;
        }
        res[i].mant = (int)res_mant_temp[i];
        int implicit_1 = (int)(pow(2, BFP_MANT_BITSIZE));
        while((res[i].mant & implicit_1) != implicit_1){   //11.01 x 2^1 = 1.101 x 2^2
            res[i].mant <<= 1;
            res[i].exp -= 1;
        }
        // cout << "\nafter normalitzation\n";
        // printBit_bfpNum(res[i], true);
    }

    // cout << "\nresult of addition\n";
    // print_bfpNum(res[0]);
    // print_bfpNum(res[1]);
    // print_bfpNum(res[2]);
    return bfpNums_to_color(res);
}


color mult_color_bfpBlock(bfpBlock block){
    vector<bfpNum> res(3, {0, (unsigned int)(block.common_exp * block.M.size() / 3  - ((int)pow(2, BFP_EXP_BITSIZE - 1) - 1) * (block.M.size()/3 - 1)), 0});
    vector<unsigned long long> res_mant_temp{(unsigned long long)block.M[0], (unsigned long long)block.M[1], (unsigned long long)block.M[2]};

    //1. assign sign
    for(int i=0; i<block.sign.size(); i+=3){
        res[0].sign ^= block.sign[i];
        res[1].sign ^= block.sign[i + 1];
        res[2].sign ^= block.sign[i + 2];
    }

    for(int i=3; i<block.M.size(); i+=3){
        //2. multiply mantissas
        res_mant_temp[0] *= (unsigned long long) block.M[i];
        res_mant_temp[1] *= (unsigned long long) block.M[i + 1];
        res_mant_temp[2] *= (unsigned long long) block.M[i + 2];

        //3 rounding
        for(int j=0; j<3; j++){
            int temp = (int)pow(2, BFP_MANT_BITSIZE); 
            unsigned short last_bit =  (unsigned short)((res_mant_temp[i] & temp) >> BFP_MANT_BITSIZE);
            temp >>= 1;
            unsigned short ground_bit = (unsigned short)((res_mant_temp[i] & temp) >> BFP_MANT_BITSIZE - 1);
            temp >>= 1;
            unsigned short round_bit =  (unsigned short)((res_mant_temp[i] & temp) >> BFP_MANT_BITSIZE - 2);
            temp -= 1;
            unsigned short sticky_bits = (unsigned short)(res_mant_temp[i] & temp);

            int lsb = (int)pow(2, BFP_MANT_BITSIZE); 
            if(ground_bit == 1){
                if(round_bit == 0 && sticky_bits == 0){ //round to even
                    if(last_bit == 1){
                        res_mant_temp[i] += lsb;
                    }
                }
                else{
                    res_mant_temp[i] += lsb; //round up
                }
            }

        }

        res_mant_temp[0] >>= BFP_MANT_BITSIZE;
        res_mant_temp[1] >>= BFP_MANT_BITSIZE;
        res_mant_temp[2] >>= BFP_MANT_BITSIZE;
    }

    //3. normalization
    for(int i=0; i<3; i++){   
        while((res_mant_temp[i] >> (BFP_MANT_BITSIZE + 1)) > 0){ //when carry is 1: 11.01 x 2^1 = 1.101 x 2^2
            res_mant_temp[i] >>= 1;
            res[i].exp += 1;
            // printf("\tCARRY! : ");
            // printBit_bfpNum_exp(res.exp, false);
            // printBit_ulong(res_mant_temp, true);
        }    
        res[i].mant = (int)res_mant_temp[i];
        int implicit_1 = (int)(pow(2, BFP_MANT_BITSIZE));
        while((res[i].mant & implicit_1) != implicit_1){   //0.101 x 2^2 = 1.01 x 2^1
            res[i].mant <<= 1;
            res[i].exp -= 1;
        }
    }

    return bfpNums_to_color(res);
}
