#pragma once
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>
#include <assert.h>

#include "bfpStruct.cuh"

using namespace std;

__managed__ unsigned int BFP_MANT_EXTRACT_BITS = 4294967295;
__managed__ unsigned int BFP_EXP_EXTRACT_BITS = 255;
__managed__ unsigned int FLOAT_MANT_EXTRACT_BITS = 0x7fffff;
__managed__ unsigned int FLOAT_EXP_EXTRACT_BITS = 0x7f800000;
__managed__ unsigned FLOAT_SIGN_EXTRACT_BITS = 0x80000000;
__managed__ unsigned long long int LONG_ALL_ONES = 0xffffffffffffffff;
__managed__ unsigned int INT_ALL_ONES = 0xffffffff;

namespace bfp
{

    /* type conversions */
    __device__ float bfpNum_to_float(bfpNum b)
    {
        unsigned int t = 0;
        int sign, exp, mant;

        sign = (int)b.sign;
        exp = (int)b.exp;
        mant = (int)b.mant;

        if (FLOAT_EXP_BITSIZE < BFP_EXP_BITSIZE)
            exp >>= (BFP_EXP_BITSIZE - FLOAT_EXP_BITSIZE);
        else
            exp <<= (FLOAT_EXP_BITSIZE - BFP_EXP_BITSIZE);

        if (FLOAT_MANT_BITSIZE < BFP_MANT_BITSIZE)
            mant >>= (BFP_MANT_BITSIZE - FLOAT_MANT_BITSIZE);
        else
            mant <<= (FLOAT_MANT_BITSIZE - BFP_MANT_BITSIZE);

        mant &= FLOAT_MANT_EXTRACT_BITS; // remove implicit 1

        t ^= sign;
        t <<= FLOAT_EXP_BITSIZE;
        t ^= exp;
        t <<= FLOAT_MANT_BITSIZE;
        t ^= mant;

        float *res = reinterpret_cast<float *>(&t);
        return *res;
    }

    __device__ bfpNum float_to_bfpNum(float f)
    {
        bfpNum res = {0};

        // cast to integer for bitwise operations
        unsigned int *t = reinterpret_cast<unsigned int *>(&f);

        // extract sign, exponent, and mantissa bits using bitwise AND operation
        res.mant = (*t & FLOAT_MANT_EXTRACT_BITS);
        res.exp = (*t & FLOAT_EXP_EXTRACT_BITS) >> FLOAT_MANT_BITSIZE;
        res.sign = (*t & FLOAT_SIGN_EXTRACT_BITS) >>
                   (FLOAT_MANT_BITSIZE + FLOAT_EXP_BITSIZE);

        // add implicit 1 when needed
        if (res.exp)
        {
            res.mant ^=
                0x00800000;
        }

        // set bfp's mantissa from float's extracted mantissa bits
        if (FLOAT_MANT_BITSIZE >=
            BFP_MANT_BITSIZE)
        { // if bfp mantissa is smaller in bit size than float
          // mantissa
            res.mant >>= (FLOAT_MANT_BITSIZE - BFP_MANT_BITSIZE);
        }

        // set bfp's exponent from float's extracted exponent bits
        if (FLOAT_EXP_BITSIZE >= BFP_EXP_BITSIZE)
        { // if bfp exponent is smaller in
          // bit size than float exponent
            res.exp >>= (FLOAT_EXP_BITSIZE - BFP_EXP_BITSIZE);
        }

        return res;
    }

    __device__ bfpNum int_to_bfpNum(int i) { return float_to_bfpNum((float)i); }

    __device__ int bfpNum_to_int(bfpNum b) { return int(bfpNum_to_float(b)); }

    /* frequently used bfpNums */
    // __managed__ bfpNum b_1 = int_to_bfpNum(1);
    __managed__ bfpNum b_0 = {0, 0, 0};

    // __managed__ bfpNum b_2 = int_to_bfpNum(2);
    // __managed__ bfpNum b_1_neg = int_to_bfpNum(-1);
    // __managed__ bfpNum b_pi = float_to_bfpNum(3.1415926535897932385);
    // __managed__ bfpNum b_180 = int_to_bfpNum(180);
    // __managed__ bfpNum b_256 = int_to_bfpNum(256);
    // __managed__ bfpNum b_infinity =
    // float_to_bfpNum(std::numeric_limits<float>::infinity());
    // __managed__ bfpNum b_0_999 = float_to_bfpNum(0.999);

    /* block formatting */
    // __device__ bfpBlock createBfpBlock(vector<float> X) {
    //     bfpBlock block;

    //     // transform floating-point numbers into bfp format
    //     vector<bfpNum> temp;
    //     for (int i = 0; i < X.size(); i++) {
    //         temp.push_back(
    //             float_to_bfpNum(X[i])); // adds implicit 1 in the process!
    //     }

    //     // find common exponent
    //     unsigned int max_exp = 0;
    //     for (int i = 0; i < X.size(); i++) {
    //         if (temp[i].exp > max_exp) {
    //             max_exp = temp[i].exp;
    //         }
    //     }
    //     block.common_exp = max_exp;

    //     // save aligned mantissas to block
    //     for (int i = 0; i < temp.size(); i++) {
    //         block.M.push_back(temp[i].mant >> (max_exp - temp[i].exp));
    //         block.sign.push_back(temp[i].sign);
    //     }

    //     return block;
    // }

    // __device__ bfpBlock createBfpBlock(vector<bfpNum> X) {
    //     bfpBlock block;

    //     // find common exponent
    //     unsigned int max_exp = 0;
    //     for (int i = 0; i < X.size(); i++) {
    //         if (X[i].exp > max_exp) {
    //             max_exp = X[i].exp;
    //         }
    //     }
    //     block.common_exp = max_exp;

    //     // save aligned mantissas to block
    //     for (int i = 0; i < X.size(); i++) {
    //         block.M.push_back(X[i].mant >> (max_exp - X[i].exp));
    //         block.sign.push_back(X[i].sign);
    //     }

    //     return block;
    // }

    /* arithmetic operations for 2 numbers */
    __device__ bfpNum add(bfpNum a, bfpNum b)
    {
        bfpNum res = {0, 0, 0};

        // if both numbers are 0, skip process
        if ((a.sign == 0) && (a.mant == 0) && (a.exp == 0) && (b.sign == 0) &&
            (b.mant == 0) && (b.exp == 0))
            return b_0;

        /* save mantissas in long long so that there is no data loss during
         * shifting */
        int temp_shift_num = 61 - BFP_MANT_BITSIZE;
        long long res_mant_temp = 0;
        long long a_mant_temp = (long long)a.mant << temp_shift_num;
        long long b_mant_temp = (long long)b.mant << temp_shift_num;
        long long implicit_1 = (long long)1 << 61;

        /* decide exponent */
        if (a.exp >= b.exp)
        {
            res.exp = a.exp;
            b_mant_temp >>= (a.exp - b.exp);
        }
        else
        {
            res.exp = b.exp;
            a_mant_temp >>= (b.exp - a.exp);
        }

        // if addition result is 0, skip process
        if ((a.sign != b.sign) && (a_mant_temp == b_mant_temp))
            return b_0;

        /* cal mantissa */
        // 1. conversion to 2's complement for negative mantissas
        if (a.sign)
        {
            a_mant_temp = ~a_mant_temp + 1;
        }
        if (b.sign)
        {
            b_mant_temp = ~b_mant_temp + 1;
        }

        // 2. add mantissas
        res_mant_temp = a_mant_temp + b_mant_temp;

        // 3. convert to signed magnitude if negative
        if (res_mant_temp & 0x8000000000000000)
        {
            res_mant_temp = ~res_mant_temp + 1;
            res.sign = 1;
        }

        // 4. normalization(implicit 1)
        if (a.sign | b.sign) // add implicit 1 for negative number addition
            res_mant_temp |= implicit_1;
        bool has_implicit_1 = (bool)(res_mant_temp >> 61);
        if (res.exp && !has_implicit_1)
        {
            while (!(res_mant_temp & implicit_1))
            { // if there is no implicit 1
                res_mant_temp <<= 1;
                res.exp -= 1;
            }
        }

        int carry = (int)(res_mant_temp >> 62);
        do
        {
            // 5. normalization(carry)
            while (carry)
            { // if there is carry
                res_mant_temp >>= 1;
                res.exp += 1;
                carry >>= 1;
            }

            // 6. rounding to nearest even
            long long t = (long long)1 << temp_shift_num;
            unsigned short last_bit =
                (unsigned short)((res_mant_temp & t) >> temp_shift_num);
            t >>= 1;
            unsigned short ground_bit =
                (unsigned short)((res_mant_temp & t) >> temp_shift_num - 1);
            t >>= 1;
            unsigned short round_bit =
                (unsigned short)((res_mant_temp & t) >> temp_shift_num - 2);
            t -= 1;
            unsigned short sticky_bits = (unsigned short)((bool)(res_mant_temp & t));

            long long lsb = (long long)1 << temp_shift_num;
            if (ground_bit)
            {
                if (round_bit == 0 && sticky_bits == 0)
                { // round to even
                    if (last_bit)
                    {
                        res_mant_temp += lsb;
                    }
                }
                else
                {
                    res_mant_temp += lsb; // round up
                }
            }

            carry = (int)(res_mant_temp >> 62);
        } while (carry);

        // 7. store result
        res.mant = (int)(res_mant_temp >> temp_shift_num);

        return res;
    }

    __device__ bfpNum sub(bfpNum a, bfpNum b)
    {
        b.sign ^= 1;
        return add(a, b);
    }

    __device__ bfpNum mult(bfpNum a, bfpNum b)
    {
        bfpNum res = {(unsigned short)(a.sign ^ b.sign), a.exp + b.exp - 127, 0};
        unsigned long long res_mant_temp = 0;

        // if a or b is 0
        if (a.mant == 0 && a.exp == 0)
            return b_0;
        if (b.mant == 0 && b.exp == 0)
            return b_0;

        // underflow: if number is too small return 0
        if ((int)res.exp < 0)
            return b_0;

        // 1. multiply
        res_mant_temp = (unsigned long long)a.mant * (unsigned long long)b.mant;

        // 3. normalization(implicit 1)
        if (res.exp)
        {
            unsigned long long implicit_1 = (unsigned long long)1 << (BFP_MANT_BITSIZE + BFP_MANT_BITSIZE);
            bool has_implicit_1 = (bool)(res_mant_temp >> (BFP_MANT_BITSIZE + BFP_MANT_BITSIZE));

            if (!has_implicit_1)
            { // exponent==0이면 implicit 1을 추가하지 않기 때문
                while (res.exp && !(res_mant_temp & implicit_1))
                { // if there is no implicit 1
                    res_mant_temp <<= 1;

                    res.exp -= 1;
                }
            }

            if (res.exp && (a.exp == 0 ^ b.exp == 0))
                res.exp += 1;
        }

        int carry = (int)(res_mant_temp >> (BFP_MANT_BITSIZE + BFP_MANT_BITSIZE + 1));
        do
        {
            // 4. normalization(Carry)
            while (carry)
            { // if there is carry
                res_mant_temp >>= 1;
                res.exp += 1;
                carry >>= 1;
            }

            // 5. rounding to nearest even
            int t = (int)std::pow(2, BFP_MANT_BITSIZE);
            unsigned short last_bit =
                (unsigned short)((res_mant_temp & t) >> BFP_MANT_BITSIZE);
            t >>= 1;
            unsigned short ground_bit =
                (unsigned short)((res_mant_temp & t) >> BFP_MANT_BITSIZE - 1);
            t >>= 1;
            unsigned short round_bit =
                (unsigned short)((res_mant_temp & t) >> BFP_MANT_BITSIZE - 2);
            t -= 1;
            unsigned short sticky_bits = (unsigned short)((bool)(res_mant_temp & t));

            int lsb = (int)std::pow(2, BFP_MANT_BITSIZE);
            if (ground_bit)
            {
                if (round_bit == 0 && sticky_bits == 0)
                { // round to even
                    if (last_bit)
                    {
                        res_mant_temp += lsb;
                    }
                }
                else
                {
                    res_mant_temp += lsb; // round up
                }
            }

            carry = (int)(res_mant_temp >> (BFP_MANT_BITSIZE + BFP_MANT_BITSIZE + 1)); // update carry
        } while (carry);

        res.mant =
            (int)(res_mant_temp >> BFP_MANT_BITSIZE); // save result in res.mant

        return res;
    }

    __device__ bfpNum div(bfpNum a, bfpNum b)
    {
        bfpNum res = {(unsigned short)(a.sign ^ b.sign), a.exp - b.exp + 127, 0};

        // if a is 0
        if (a.sign == 0 && a.exp == 0 && a.mant == 0)
            return b_0;

        // underflow: if number is too small return 0
        if ((int)res.exp < 0)
            return b_0;

        // 0. conversion to temps & shifting
        unsigned long long a_temp = (unsigned long long)a.mant;
        unsigned long long b_temp = (unsigned long long)b.mant;
        unsigned long long msb = (unsigned long long)1 << 63;
        long long shift_cnt = -(64 - BFP_SIGNIFICAND_BITSIZE);

        while (!(a_temp & msb))
        {
            shift_cnt++;
            a_temp <<= 1;
        }
        res.exp = max((long long)0, res.exp - shift_cnt);

        // 1. divide mantissas
        unsigned long long res_mant_temp = (unsigned long long)a_temp / b_temp;

        // 2. normalization(implicit 1)
        if (res.exp)
        {
            unsigned long long implicit_1 = (unsigned long long)1 << (64 - BFP_SIGNIFICAND_BITSIZE);
            bool has_implicit_1 = (bool)(res_mant_temp >> (64 - BFP_SIGNIFICAND_BITSIZE));

            if (!has_implicit_1)
            {
                while (!(res_mant_temp & implicit_1))
                {
                    res_mant_temp <<= 1;
                    res.exp -= 1;
                }
            }

            if (a.exp == 0 ^ b.exp == 0)
                res.exp += 1;
        }
        else
        {
            res_mant_temp >>= 1;
            // res_mant_temp += (unsigned long long)1 << (64 - 1 - BFP_MANT_BITSIZE - BFP_MANT_BITSIZE);
        }

        int carry = (int)(res_mant_temp >> (64 - BFP_SIGNIFICAND_BITSIZE + 1));
        do
        {
            // 3. normalization(carry)
            while (carry)
            { // if there is carry
                res_mant_temp >>= 1;
                res.exp += 1;
                carry >>= 1;
            }

            // 4. rounding to nearest even
            int lsb_zero_cnt = 64 - 1 - BFP_MANT_BITSIZE - BFP_MANT_BITSIZE;
            unsigned long long t = (unsigned long long)1 << lsb_zero_cnt;
            unsigned short last_bit =
                (unsigned short)((res_mant_temp & t) >> lsb_zero_cnt);
            t >>= 1;
            unsigned short ground_bit =
                (unsigned short)((res_mant_temp & t) >> (lsb_zero_cnt - 1));
            t >>= 1;
            unsigned short round_bit =
                (unsigned short)((res_mant_temp & t) >> (lsb_zero_cnt - 2));
            t -= 1;
            unsigned short sticky_bits = (unsigned short)((bool)(res_mant_temp & t));

            unsigned long long lsb = (unsigned long long)1 << lsb_zero_cnt;
            if (ground_bit)
            {
                if (round_bit == 0 && sticky_bits == 0)
                {
                    if (last_bit)
                    {
                        res_mant_temp += lsb;
                    }
                }
                else
                {
                    res_mant_temp += lsb;
                }
            }

            carry = (int)(res_mant_temp >> (64 - BFP_SIGNIFICAND_BITSIZE + 1));

        } while (carry);

        int lsb_zero_cnt = 64 - 1 - BFP_MANT_BITSIZE - BFP_MANT_BITSIZE;
        res.mant = (int)(res_mant_temp >> lsb_zero_cnt); // save result in res.mant

        return res;
    }

    __device__ bfpNum sqrt(bfpNum a)
    {
        return float_to_bfpNum(sqrtf(bfpNum_to_float(a)));
    }

    __device__ bfpNum sqrt(float a) {return float_to_bfpNum(sqrtf(a)); }

    __device__ bfpNum pow(bfpNum base, bfpNum n)
    {
        return float_to_bfpNum(powf(bfpNum_to_float(base), bfpNum_to_float(n)));
    }

    __device__ bfpNum pow(bfpNum base, float n)
    {
        return float_to_bfpNum(powf(bfpNum_to_float(base), n));
    }

    __device__ bfpNum abs(bfpNum a)
    {
        a.sign = 0;
        return a;
    }

    __device__ bfpNum tan(bfpNum a)
    {
        return float_to_bfpNum(tanf(bfpNum_to_float(a)));
    }

    __device__ bfpNum sin(bfpNum a)
    {
        return float_to_bfpNum(sinf(bfpNum_to_float(a)));
    }

    __device__ bfpNum cos(bfpNum a)
    {
        return float_to_bfpNum(cosf(bfpNum_to_float(a)));
    }

    __device__ bool compare(bfpNum a, bfpNum b)
    {
        return bfpNum_to_float(a) > bfpNum_to_float(b);
    }

    __device__ bool isequal(bfpNum a, bfpNum b)
    {
        return bfpNum_to_float(a) == bfpNum_to_float(b);
    }

    __device__ inline bfpNum operator+(bfpNum a, bfpNum b) { return add(a, b); }

    __device__ inline bfpNum operator-(bfpNum a, bfpNum b) { return sub(a, b); }

    __device__ inline bfpNum operator-(bfpNum a)
    {
        a.sign ^= 1;
        return a;
    }

    __device__ inline bfpNum operator*(bfpNum a, bfpNum b) { return mult(a, b); }

    __device__ inline bfpNum operator/(bfpNum a, bfpNum b) { return div(a, b); }

    __device__ inline bfpNum min(bfpNum a, bfpNum b)
    {
        float a_f = bfpNum_to_float(a);
        float b_f = bfpNum_to_float(b);
        float min = fminf(a_f, b_f);
        return float_to_bfpNum(min);
    }

    __device__ inline bfpNum max(bfpNum a, bfpNum b)
    {
        float a_f = bfpNum_to_float(a);
        float b_f = bfpNum_to_float(b);
        float max = fmaxf(a_f, b_f);
        return float_to_bfpNum(max);
    }
    
    __device__ inline bool operator>(bfpNum a, bfpNum b) { return compare(a, b); }

    __device__ inline bool operator>=(bfpNum a, bfpNum b)
    {
        return compare(a, b) || isequal(a, b);
    }

    __device__ inline bool operator<(bfpNum a, bfpNum b) { return compare(b, a); }

    __device__ inline bool operator<=(bfpNum a, bfpNum b)
    {
        return compare(b, a) || isequal(a, b);
    }

    __device__ inline bool operator==(bfpNum a, bfpNum b) { return isequal(a, b); }

    __device__ inline bool operator!=(bfpNum a, bfpNum b)
    {
        return (a.sign != b.sign) || (a.exp != b.exp) || (a.mant != b.mant);
    }

    /* arithmetic operations for entire block */
    // __device__ bfpNum add_bfpBlock_b(bfpBlock block) {
    //     bfpNum res = {0, block.common_exp, 0};
    //     long long int res_mant_temp = 0;

    //     // 1. converision to 2's complment for negative mantissas
    //     for (int i = 0; i < block.sign.size(); i++) {
    //         if (block.sign[i]) {
    //             block.M[i] = ~block.M[i] + 1;
    //         }
    //     }

    //     // 2. add mantissas
    //     for (int i = 0; i < block.M.size(); i++) {
    //         res_mant_temp += (long long int)block.M[i];
    //         // printf("add ");
    //         // printBit_sint(block.M[i], false);
    //         // printf(" => ");
    //         // printBit_sint(res.mant, true);

    //         /* option 1 */
    //         // if((res.mant & 0xff000000) > 0){ //when carry is 1: 11.01 x 2^1
    //         // = 1.101 x 2^2
    //         //     res.mant >>= 1;
    //         //     res.exp += 1;
    //         //     printf("\tCARRY! : ");
    //         //     printBit_mant(res.mant, true);
    //         // }
    //     }

    //     // 3. convert to  signed magnitude if negative
    //     if ((res_mant_temp & 0x8000000000000000) ==
    //         0x8000000000000000) { // if MSB is 1(sign = 1)
    //         res_mant_temp = ~res_mant_temp + 1;
    //         res.sign = 1;
    //     }

    //     // // if addition result is 0
    //     // if(res_mant_temp == 0){

    //     // }

    //     // 4. normalization
    //     /* option 2 */
    //     while ((res_mant_temp >> (BFP_MANT_BITSIZE + 1)) >
    //            0) { // when carry is 1: 11.01 x 2^1 = 1.101 x 2^2
    //         res_mant_temp >>= 1;
    //         res.exp += 1;
    //         // printf("\tCARRY! : ");
    //         // printBit_bfpNum_exp(res.exp, false);
    //         // printBit_ulong(res_mant_temp, true);
    //     }
    //     res.mant = (int)res_mant_temp;
    //     int implicit_1 = (int)(std::pow(2, BFP_MANT_BITSIZE));
    //     while ((res.mant & implicit_1) != implicit_1) { // 0.101 x 2^2 = 1.01 x
    //     2^1
    //         res.mant <<= 1;
    //         res.exp -= 1;
    //     }

    //     return res;
    // }

    // __device__ float add_bfpBlock(bfpBlock block) {
    //     return bfpNum_to_float(add_bfpBlock_b(block));
    // }

    // __device__ bfpNum mult_bfpBlock_b(bfpBlock block) {
    //     bfpNum res = {0,
    //                   (unsigned int)(block.common_exp * block.M.size() -
    //                                  ((int)std::pow(2, BFP_EXP_BITSIZE - 1) - 1)
    //                                  *
    //                                      (block.M.size() - 1)),
    //                   0};

    //     // 1. assign sign
    //     for (int i = 0; i < block.sign.size(); i++) {
    //         res.sign ^= block.sign[i];
    //     }

    //     unsigned long long res_mant_temp = (unsigned long long)block.M[0];

    //     for (int i = 1; i < block.M.size(); i++) {
    //         // 2. multiply mantissas
    //         res_mant_temp *= (unsigned long long)block.M[i];

    //         // 3 rounding
    //         int temp = (int)std::pow(2, BFP_MANT_BITSIZE);
    //         unsigned short last_bit =
    //             (unsigned short)((res_mant_temp & temp) >> BFP_MANT_BITSIZE);
    //         temp >>= 1;
    //         unsigned short ground_bit =
    //             (unsigned short)((res_mant_temp & temp) >> BFP_MANT_BITSIZE - 1);
    //         temp >>= 1;
    //         unsigned short round_bit =
    //             (unsigned short)((res_mant_temp & temp) >> BFP_MANT_BITSIZE - 2);
    //         temp -= 1;
    //         unsigned short sticky_bits = (unsigned short)(res_mant_temp & temp);

    //         int lsb = (int)std::pow(2, BFP_MANT_BITSIZE);
    //         if (ground_bit == 1) {
    //             if (round_bit == 0 && sticky_bits == 0) { // round to even
    //                 if (last_bit == 1) {
    //                     res_mant_temp += lsb;
    //                 }
    //             } else {
    //                 res_mant_temp += lsb; // round up
    //             }
    //         }

    //         res_mant_temp >>= BFP_MANT_BITSIZE;
    //     }

    //     // 3. normalization
    //     while ((res_mant_temp >> (BFP_MANT_BITSIZE + 1)) >
    //            0) { // when carry is 1: 11.01 x 2^1 = 1.101 x 2^2
    //         res_mant_temp >>= 1;
    //         res.exp += 1;
    //         // printf("\tCARRY! : ");
    //         // printBit_bfpNum_exp(res.exp, false);
    //         // printBit_ulong(res_mant_temp, true);
    //     }
    //     res.mant = (int)res_mant_temp;
    //     int implicit_1 = (int)(std::pow(2, BFP_MANT_BITSIZE));
    //     while ((res.mant & implicit_1) != implicit_1) { // 0.101 x 2^2 = 1.01 x
    //     2^1
    //         res.mant <<= 1;
    //         res.exp -= 1;
    //     }

    //     return res;
    // }

    // __device__ float mult_bfpBlock(bfpBlock block) {
    //     return bfpNum_to_float(mult_bfpBlock_b(block));
    // }

} // namespace bfp
