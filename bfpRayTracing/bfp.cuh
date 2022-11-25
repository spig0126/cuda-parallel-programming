#include <float.h>

#include "bfpStruct.cuh"
#include "print.cuh"
#include "vec3.h"

using namespace std;


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

vector<bfpNumFloat> color_to_bfpNumFloat(color c){
    vector<bfpNumFloat> c_formatted;

    for(int i=0; i<3; i++){
        c_formatted.push_back(float_to_bfpNumFloat(c[i]));
    }

    return c_formatted;
}

color bfpNumFloats_to_color(vector<bfpNumFloat> v){
    float r = bfpNumFloat_to_float(v[0]);
    float g = bfpNumFloat_to_float(v[1]);
    float b = bfpNumFloat_to_float(v[2]);

    return color(r, g, b);
}



/* block formatting */
bfpBlock createBfpBlock(vector<float> X){
    bfpBlock block;

    //transform floating-point numbers into bfp format
    vector<bfpNumFloat> temp;
    for(int i=0; i<X.size(); i++){
        temp.push_back(float_to_bfpNumFloat(X[i])); //adds implicit 1 in the process!
    } 

    //find common exponent
    unsigned int max_exp = 0;
    for(int i=0; i<X.size(); i++){
        if(temp[i].exp > max_exp){
            max_exp = temp[i].exp;
        }
    }

    //save common exponent
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
    vector<bfpNumFloat> temp;
    for(int i=0; i<colors.size(); i++){
        vector<bfpNumFloat> c_formatted = color_to_bfpNumFloat(colors[i]);  //adds implicit 1 in the process!
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
bfpNumFloat div_f(bfpBlock block, bfpNumFloat a, bfpNumFloat b){
    bfpNumFloat res = {(unsigned short)(a.sign^b.sign), 127, 0};

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
    bfpNumFloat res = {0, block.common_exp, 0};

    //1. converision to 2's complment for negative mantissas
    for(int i=0; i<block.sign.size(); i++){
        if(block.sign[i]){
            block.M[i] = ~block.M[i] + 1;
        }
    }

    //2. add mantissas
    for(int i=0; i<block.M.size(); i++){       
        res.mant += block.M[i];
    }

    //3. convert to  signed magnitude if negative
    if((res.mant & 0x80000000) == 0x80000000){    //if MSB is 1(sign = 1)
        res.mant = ~res.mant + 1;
        res.sign = 1;
    }

    //4. normalization
    while(res.mant & 0x01000000){ //when carry is 1
        res.mant >>= 1;
        res.exp += 1;
    }
    while((res.mant & 0x00800000) != 0x00800000){   //11.01 x 2^1 = 1.101 x 2^2
        res.mant <<= 1;
        res.exp -= 1;
    }

    res.mant &= 0x007fffff;

    return bfpNumFloat_to_float(res);
}

float mult_bfpBlock(bfpBlock block){
    bfpNumFloat res = {0, (unsigned int)(block.common_exp * block.M.size() - 127 * (block.M.size()- 1)), 0};

    //1. assign sign
    for(int i=0; i<block.sign.size(); i++){
        res.sign ^= block.sign[i];
    }

    unsigned long long res_temp = (unsigned long long)block.M[0];

    for(int i=1; i<block.M.size(); i++){
        //2. multiply mantissas
        res_temp *= (unsigned long long) block.M[i];

        //3 rounding
        unsigned short last_bit =  (unsigned short)((res_temp & 0x0000000000800000) >> 23);
        unsigned short ground_bit = (unsigned short)((res_temp & 0x0000000000400000) >> 22);
        unsigned short round_bit =  (unsigned short)((res_temp & 0x0000000000200000) >> 21);
        unsigned short sticky_bits = (unsigned short)(res_temp & 0x00000000001fffff);

        if(ground_bit == 1){
            if(round_bit == 0 && sticky_bits == 0){ //round to even
                if(last_bit == 1){
                    res_temp += 0x000000000800000;
                }
            }
            else{
                res_temp += 0x000000000800000; //round up
            }
        }

        res_temp >>= 23;
    }


    // 3. normalization
    while((res_temp & 0x800000) != 0x800000 && (res_temp != 0)){
        res_temp <<= 1;
        res.exp-=1;
    }  
    // int temp = 0;
    // while((res_temp & 0xffffffffff000000) > 0){
    //     res_temp >>= 1;
    //     temp++;
    // } 
    // res.exp-= (temp - 23) > 0? temp - 23 : 0;


    //5. normalization(carry)
    int carry = (int)res_temp >> 24;
    while(carry > 0){
        res_temp >>= 1;
        carry >>= 1;
        res.exp += 1;
    }

    //6. remove implicit 1
    res.mant = (int)res_temp & 0x007fffff;

    
    return bfpNumFloat_to_float(res);
}

//color block
color add_color_bfpBlock(bfpBlock block){
    vector<bfpNumFloat> res(3, {0,block.common_exp,0});

    //1. converision to 2's complment for negative mantissas
    for(int i=0; i<block.sign.size(); i++){
        if(block.sign[i]){
            block.M[i] = ~block.M[i] + 1;
        }
    }

    //2. add mantissas
    for(int i=0; i<block.M.size(); i+=3){
        res[0].mant += block.M[i];
        res[1].mant += block.M[i + 1];
        res[2].mant += block.M[i + 2];

        // printBit_sint(res[0].mant, true);
        // printBit_sint(res[1].mant, true);
        // printBit_sint(res[2].mant, true);
    }

    //3. convert to  signed magnitude if negative
    for(int i=0; i<3; i++){
        if((res[i].mant & 0x80000000) == 0x80000000){    //if MSB is 1(sign = 1)
            res[i].mant = ~res[i].mant + 1;
            res[i].sign = 1;
        }
    }

    //4. normalization
    for(int i=0; i<3; i++){
        while(res[i].mant & 0x01000000){ //when carry is 1
            res[i].mant >>= 1;
            res[i].exp += 1;
        }
        while((res[i].mant & 0x00800000) != 0x00800000){   //11.01 x 2^1 = 1.101 x 2^2
            res[i].mant <<= 1;
            res[i].exp -= 1;
        }
        // cout << "\nafter normalitzation\n";
        // printBit_bfpNumFloat(res[i], true);
    }

    //5. remove implicit 1
    for(int i=0; i<3; i++){
        res[i].mant &= 0x007fffff;
    }

    // cout << "\nresult of addition\n";
    // print_bfpNumFloat(res[0]);
    // print_bfpNumFloat(res[1]);
    // print_bfpNumFloat(res[2]);
    return bfpNumFloats_to_color(res);
}


color mult_color_bfpBlock(bfpBlock block){
    vector<bfpNumFloat> res(3, {0, (unsigned int)(block.common_exp * block.M.size() / 3  - 127 * (block.M.size()/3 - 1)), 0});

    //1. assign sign
    for(int i=0; i<block.sign.size(); i+=3){
        res[0].sign ^= block.sign[i];
        res[1].sign ^= block.sign[i + 1];
        res[2].sign ^= block.sign[i + 2];
    }

    vector<unsigned long long> res_temp{(unsigned long long)block.M[0], (unsigned long long)block.M[1], (unsigned long long)block.M[2]};

    for(int i=3; i<block.M.size(); i+=3){
        //2. multiply mantissas
        res_temp[0] *= (unsigned long long) block.M[i];
        res_temp[1] *= (unsigned long long) block.M[i + 1];
        res_temp[2] *= (unsigned long long) block.M[i + 2];

        //3 rounding
        for(int j=0; j<3; j++){
            unsigned short last_bit =  (unsigned short)((res_temp[j] & 0x0000000000800000) >> 22);
            unsigned short ground_bit = (unsigned short)((res_temp[j] & 0x0000000000400000) >> 22);
            unsigned short round_bit =  (unsigned short)((res_temp[j] & 0x0000000000200000) >> 21);
            unsigned short sticky_bits = (unsigned short)(res_temp[j] & 0x00000000001fffff);

            if(ground_bit == 1){
                if(round_bit == 0 && sticky_bits == 0){ //round to even
                    if(last_bit == 1){
                        res_temp[j] += 0x000000000800000;
                    }
                }
                else{
                    res_temp[j] += 0x000000000800000; //round up
                }
            }
        }

        res_temp[0] >>= 23;
        res_temp[1] >>= 23;
        res_temp[2] >>= 23;
    }


    //3. normalization
    for(int i=0; i<3; i++){
        while((res_temp[i] & 0x800000) != 0x800000 && (res_temp[i] != 0)){
            res_temp[i] <<= 1;
            res[i].exp-=1;
        }  
    }

    //5. normalization(carry)
    for(int i=0; i<3; i++){
        int carry = (int)res_temp[i] >> 24;
        while(carry > 0){
            res[i].mant >>= 1;
            carry >>= 1;
            res[i].exp += 1;
        }
    }

    //6. remove implicit 1
    for(int i=0; i<3; i++){
        res[i].mant = (int)res_temp[i] & 0x007fffff;
    }

     return bfpNumFloats_to_color(res);
}