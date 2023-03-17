#pragma once
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

#include "bfpStruct.cuh"

using namespace std;

namespace bfp
{
    int BFP_MANT_EXTRACT_BITS = (int)std::pow(2, BFP_MANT_BITSIZE) - 1;
    int BFP_EXP_EXTRACT_BITS = (int)std::pow(2, BFP_EXP_BITSIZE) - 1;
    int BFP_SIGN_EXTRACT_BITS = 1 << (BFP_MANT_BITSIZE + BFP_EXP_BITSIZE);
    int FLOAT_MANT_EXTRACT_BITS = 0x7fffff;
    int FLOAT_EXP_EXTRACT_BITS = 0x7f800000;
    int FLOAT_SIGN_EXTRACT_BITS = 0x80000000;
    long long int LONG_ALL_ONES = 0xffffffffffffffff;
    int INT_ALL_ONES = 0xffffffff;

    //---------------------------------------------------------
    void printBit_ulong(long long num, bool newLine)
    { // 4bits grouped together
        long long temp = num;
        char out[MAX_BITSIZE] = "";
        for (int i = 0; i < 64; i++, temp >>= 1)
        {
            if (i % 4 == 0)
            {
                strcat(out, " ");
            }
            if (temp & 1)
            {
                strcat(out, "1");
            }
            else
            {
                strcat(out, "0");
            }
        }

        for (int i = 79; i >= 0; i--)
        {
            std::cout << out[i];
        }

        printf("\t(%llx)", num);

        if (newLine)
        {
            std::cout << endl;
        }
    }

    void printBit_uint(unsigned int num, int len)
    {
        char out[MAX_BITSIZE] = "";

        for (int i = 0; i < len; i++, num >>= 1)
        {
            if (num & 1)
            {
                strcat(out, "1");
            }
            else
            {
                strcat(out, "0");
            }
        }

        for (int i = len - 1; i >= 0; i--)
        {
            cout << out[i];
        }
    }
    void printBit_bfpNum(bfpNum b, bool nextLine)
    {
        printBit_uint(b.sign, 1);
        cout << " ";
        printBit_uint(b.exp, BFP_EXP_BITSIZE);
        cout << " ";
        printBit_uint(b.mant, BFP_MANT_BITSIZE);

        if (nextLine)
        {
            cout << endl;
        }
    }
    //---------------------------------------------------------

    /* type conversions */
    float bfpNum_to_float(bfpNum b)
    {
        unsigned int t = 0;
        int sign, exp, mant;

        sign = (int)b.sign;
        exp = (int)b.exp;
        mant = (int)b.mant;

        if (FLOAT_EXP_BITSIZE < BFP_EXP_BITSIZE)
            exp >>= std::abs(BFP_EXP_BITSIZE - FLOAT_EXP_BITSIZE);
        else
            exp <<= std::abs(FLOAT_EXP_BITSIZE - BFP_EXP_BITSIZE);

        if (FLOAT_MANT_BITSIZE < BFP_MANT_BITSIZE)
            mant >>= std::abs(BFP_MANT_BITSIZE - FLOAT_MANT_BITSIZE);
        else
            mant <<= std::abs(FLOAT_MANT_BITSIZE - BFP_MANT_BITSIZE);

        mant &= FLOAT_MANT_EXTRACT_BITS; // remove implicit 1

        t ^= sign;
        t <<= FLOAT_EXP_BITSIZE;
        t ^= exp;
        t <<= FLOAT_MANT_BITSIZE;
        t ^= mant;

        float *res = reinterpret_cast<float *>(&t);
        return *res;
    }

    bfpNum float_to_bfpNum(float f)
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
                1 << FLOAT_MANT_BITSIZE;
        }

        // set bfp's mantissa from float's extracted mantissa bits
        if (FLOAT_MANT_BITSIZE >
            BFP_MANT_BITSIZE)
        { // if bfp mantissa is smaller in bit size than float
          // mantissa
            res.mant >>= (FLOAT_MANT_BITSIZE - BFP_MANT_BITSIZE);
        }

        // set bfp's exponent from float's extracted exponent bits
        if (FLOAT_EXP_BITSIZE > BFP_EXP_BITSIZE)
        { // if bfp exponent is smaller in
          // bit size than float exponent
            res.exp >>= (FLOAT_EXP_BITSIZE - BFP_EXP_BITSIZE);
        }


        return res;
    }

    bfpNum int_to_bfpNum(int i) { return float_to_bfpNum((float)i); }

    int bfpNum_to_int(bfpNum b) { return int(bfpNum_to_float(b)); }

    /* frequently used bfpNums */
    bfpNum b_pi = float_to_bfpNum(3.1415926535897932385);
    bfpNum b_infinity = float_to_bfpNum(std::numeric_limits<float>::infinity());

    bfpNum b_0 = {0, 0, 0};
    bfpNum b_0_001 = float_to_bfpNum(0.001);
    bfpNum b_0_1 = float_to_bfpNum(0.1);
    bfpNum b_0_15 = float_to_bfpNum(0.15);
    bfpNum b_0_2 = float_to_bfpNum(0.2);
    bfpNum b_0_3 = float_to_bfpNum(0.3);
    bfpNum b_0_4 = float_to_bfpNum(0.4);
    bfpNum b_0_5 = float_to_bfpNum(0.5);
    bfpNum b_0_6 = float_to_bfpNum(0.6);
    bfpNum b_0_7 = float_to_bfpNum(0.7);
    bfpNum b_0_8 = float_to_bfpNum(0.8);
    bfpNum b_0_9 = float_to_bfpNum(0.9);
    bfpNum b_0_95 = float_to_bfpNum(0.95);
    bfpNum b_0_999 = float_to_bfpNum(0.999);
    bfpNum b_1 = int_to_bfpNum(1);
    bfpNum b_1_neg = int_to_bfpNum(-1);
    bfpNum b_1_5 = float_to_bfpNum(1.5);
    bfpNum b_2 = int_to_bfpNum(2);
    bfpNum b_3 = int_to_bfpNum(3);
    bfpNum b_4 = int_to_bfpNum(4);
    bfpNum b_9 = int_to_bfpNum(9);
    bfpNum b_10 = int_to_bfpNum(10);
    bfpNum b_11 = int_to_bfpNum(11);
    bfpNum b_13 = int_to_bfpNum(13);
    bfpNum b_16 = int_to_bfpNum(16);
    bfpNum b_20 = int_to_bfpNum(20);
    bfpNum b_121 = int_to_bfpNum(121);
    bfpNum b_180 = int_to_bfpNum(180);
    bfpNum b_256 = int_to_bfpNum(256);
    bfpNum b_1000 = int_to_bfpNum(1000);
    bfpNum b_1000_neg = int_to_bfpNum(-1000);

    /* block formatting */
    bfpBlock createBfpBlock(vector<float> X)
    {
        bfpBlock block;

        // transform floating-point numbers into bfp format
        vector<bfpNum> temp;
        for (int i = 0; i < X.size(); i++)
        {
            temp.push_back(
                float_to_bfpNum(X[i])); // adds implicit 1 in the process!
        }

        // find common exponent
        unsigned int max_exp = 0;
        for (int i = 0; i < X.size(); i++)
        {
            if (temp[i].exp > max_exp)
            {
                max_exp = temp[i].exp;
            }
        }
        block.common_exp = max_exp;

        // save aligned mantissas to block
        for (int i = 0; i < temp.size(); i++)
        {
            block.M.push_back(temp[i].mant >> (max_exp - temp[i].exp));
            block.sign.push_back(temp[i].sign);
        }

        return block;
    }

    bfpBlock createBfpBlock(vector<bfpNum> X)
    {
        bfpBlock block;

        // find common exponent
        unsigned int max_exp = 0;
        for (int i = 0; i < X.size(); i++)
        {
            if (X[i].exp > max_exp)
            {
                max_exp = X[i].exp;
            }
        }
        block.common_exp = max_exp;

        // save aligned mantissas to block
        for (int i = 0; i < X.size(); i++)
        {
            block.M.push_back(X[i].mant >> (max_exp - X[i].exp));
            block.sign.push_back(X[i].sign);
        }

        return block;
    }

    /* arithmetic operations for 2 numbers */
    bfpNum add(bfpNum a, bfpNum b)
    {
        // std::cout << "add " << bfpNum_to_float(a) << " " << bfpNum_to_float(b) << endl;

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

        // std::cout << "conversion" << endl;
        // printBit_ulong(a_mant_temp, true);
        // printBit_ulong(b_mant_temp, true);

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

        // std::cout << "\nstep 0: decide exponent" << endl;
        // printBit_ulong(a_mant_temp, true);
        // printBit_ulong(b_mant_temp, true);
        // printBit_ulong(res.exp, true);

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

        // std::cout << "\nstep 1: conversion to 2's complement for negative mantissas"
        //           << endl;
        // printBit_ulong(a_mant_temp, true);
        // printBit_ulong(b_mant_temp, true);

        // 2. add mantissas
        res_mant_temp = a_mant_temp + b_mant_temp;

        // std::cout << "\nstep 2: add mantissas" << endl;
        // printBit_ulong(res_mant_temp, true);

        // 3. convert to signed magnitude if negative
        if (res_mant_temp & 0x8000000000000000)
        {
            res_mant_temp = ~res_mant_temp + 1;
            res.sign = 1;
        }

        // std::cout << "\nstep 3: convert to signed magnitude if negative" << endl;
        // printBit_ulong(res_mant_temp, true);

        // 4. normalization(implicit 1)
        // if (res.sign && (a.sign ^ b.sign)) // add implicit 1 for negative number addition
        //     res_mant_temp |= implicit_1;
        bool has_implicit_1 = (bool)(res_mant_temp >> 61);
        if (res.exp && !has_implicit_1)
        {
            while (!(res_mant_temp & implicit_1))
            { // if there is no implicit 1
                res_mant_temp <<= 1;
                res.exp -= 1;
            }
        }

        // std::cout << "\nstep 5: normalization(implicit 1)" << endl;
        // printBit_ulong(res_mant_temp, true);

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

            // std::cout << "\nstep 4: normalization(carry)" << endl;
            // printBit_ulong(res_mant_temp, true);

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

            // std::cout << "\nstep 6: rounding" << endl;
            // printBit_ulong(res_mant_temp, true);
            // std::cout << last_bit << ground_bit << round_bit << sticky_bits << endl;

            carry = (int)(res_mant_temp >> 62);
        } while (carry);

        // 7. store result
        res.mant = (int)(res_mant_temp >> temp_shift_num);

        // std::cout << "done" << endl;
        return res;
    }

    bfpNum sub(bfpNum a, bfpNum b)
    {
        b.sign ^= 1;
        return add(a, b);
    }

    bfpNum mult(bfpNum a, bfpNum b)
    {
        bfpNum res = {(unsigned short)(a.sign ^ b.sign), a.exp + b.exp - BFP_BIAS, 0};
        unsigned long long res_mant_temp = 0;

        // if a or b is 0
        if (a.mant == 0 && a.exp == 0)
            return b_0;
        if (b.mant == 0 && b.exp == 0)
            return b_0;

        // underflow: if number is too small return 0
        if ((int)res.exp < 0)
            return b_0;

        // std::cout << "\nstep 0: conversion to temps"
        //           << endl;
        // printBit_ulong((unsigned long long)a.mant, true);
        // printBit_ulong((unsigned long long)b.mant, true);
        // printBit_ulong((unsigned long long)res.exp, true);

        // 1. multiply
        res_mant_temp = (unsigned long long)a.mant * (unsigned long long)b.mant;

        // std::cout << "\nstep 1: multiply" << endl;
        // printBit_ulong(res_mant_temp, true);
        // printBit_ulong((unsigned long long)res.exp, true);

        // 3. normalization(implicit 1)
        if (res.exp)
        {
            unsigned long long implicit_1 = (unsigned long long)1 << (BFP_MANT_BITSIZE + BFP_MANT_BITSIZE);
            bool has_implicit_1 = (bool)(res_mant_temp >> (BFP_MANT_BITSIZE + BFP_MANT_BITSIZE));
            // long long left_shift_cnt = 0;

            if (!has_implicit_1)
            { // exponent==0이면 implicit 1을 추가하지 않기 때문
                while (res.exp && !(res_mant_temp & implicit_1))
                { // if there is no implicit 1
                    res_mant_temp <<= 1;
                    // left_shift_cnt++;
                    res.exp -= 1;
                }
            }

            if (res.exp && (a.exp == 0 ^ b.exp == 0))
                res.exp += 1;

            // res.exp = (unsigned int)(std::max((long long)0, (long long)res.exp - left_shift_cnt));
        }

        // std::cout << "\nstep 3: normalization(implicit 1)" << endl;
        // printBit_ulong(res_mant_temp, true);
        // printBit_ulong((unsigned long long)res.exp, true);

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

            // std::cout << "\nstep 4: normalization(carry)" << endl;
            // printBit_ulong(res_mant_temp, true);
            // printBit_ulong((unsigned long long)res.exp, true);

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

            // std::cout << "\nstep 5: rounding" << endl;
            // printBit_ulong(res_mant_temp, true);
            // std::cout << last_bit << ground_bit << round_bit << sticky_bits << endl;

            carry = (int)(res_mant_temp >> (BFP_MANT_BITSIZE + BFP_MANT_BITSIZE + 1)); // update carry
        } while (carry);

        res.mant =
            (int)(res_mant_temp >> BFP_MANT_BITSIZE); // save result in res.mant

        return res;
    }

    bfpNum div(bfpNum a, bfpNum b)
    {

        bfpNum res = {(unsigned short)(a.sign ^ b.sign), a.exp - b.exp + BFP_BIAS, 0};

        // if a is 0
        if (a.sign == 0 && a.exp == 0 && a.mant == 0)
            return b_0;
        if (b.sign == 0 && b.exp == 0 && b.mant == 0)
        {
            return b_infinity;
            // throw std::invalid_argument("EXCEPTION: division with 0");
            // exit(0);
        }

        // underflow: if number is too small return 0
        if ((int)res.exp < 0)
            return b_0;

        // 0. conversion to temps & shifting
        unsigned long long a_temp = (unsigned long long)a.mant;
        unsigned long long b_temp = (unsigned long long)b.mant;

        // std::cout << "\nstep 0: conversion to temps"
        //           << endl;
        // printBit_ulong(a_temp, true);
        // printBit_ulong(b_temp, true);
        // printBit_ulong((unsigned long long)res.exp, true);

        unsigned long long msb = (unsigned long long)1 << 63;
        long long shift_cnt = -(64 - BFP_SIGNIFICAND_BITSIZE);

        while (!(a_temp & msb))
        {
            shift_cnt++;
            a_temp <<= 1;
        }
        res.exp = std::max((long long)0, res.exp - shift_cnt);

        // // unsigned long long a_temp = (unsigned long long)a.mant << (64 - BFP_SIGNIFICAND_BITSIZE);

        // cout << "shift cnt: " << shift_cnt << endl;
        // printBit_ulong(a_temp, true);
        // printBit_ulong(b_temp, true);
        // printBit_ulong((unsigned long long)res.exp, true);

        // 1. divide mantissas
        unsigned long long res_mant_temp = (unsigned long long)a_temp / b_temp;

        // std::cout << "\nstep 1: division" << endl;
        // printBit_ulong(res_mant_temp, true);
        // printBit_ulong((unsigned long long)res.exp, true);

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

        // std::cout << "\nstep 2: normalization(implicit 1)" << endl;
        // printBit_ulong(res_mant_temp, true);
        // printBit_ulong((unsigned long long)res.exp, true);

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

            // std::cout << "\nstep 4: normalization(carry)" << endl;
            // printBit_ulong(res_mant_temp, true);
            // printBit_ulong((unsigned long long)res.exp, true);

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

            // std::cout << "\nstep 4: rounding" << endl;
            // printBit_ulong(res_mant_temp, true);
            // printBit_ulong((unsigned long long)res.exp, true);
            // std::cout << last_bit << ground_bit << round_bit << sticky_bits << endl;

            carry = (int)(res_mant_temp >> (64 - BFP_SIGNIFICAND_BITSIZE + 1)); 

        } while (carry);

        int lsb_zero_cnt = 64 - 1 - BFP_MANT_BITSIZE - BFP_MANT_BITSIZE;
        res.mant = (int)(res_mant_temp >> lsb_zero_cnt); // save result in res.mant

        return res;
    }

    bfpNum sqrt(bfpNum a) { return float_to_bfpNum(std::sqrt(bfpNum_to_float(a))); }

    bfpNum sqrt(float a) { return float_to_bfpNum(std::sqrt(a)); }

    bfpNum pow(bfpNum base, bfpNum n)
    {
        return float_to_bfpNum(std::pow(bfpNum_to_float(base), bfpNum_to_float(n)));
    }

    bfpNum pow(bfpNum base, float n)
    {
        return float_to_bfpNum(std::pow(bfpNum_to_float(base), n));
    }

    bfpNum abs(bfpNum a)
    {
        a.sign = 0;
        return a;
    }

    bfpNum tan(bfpNum a) { return float_to_bfpNum(std::tan(bfpNum_to_float(a))); }

    bfpNum sin(bfpNum a) { return float_to_bfpNum(std::sin(bfpNum_to_float(a))); }

    bfpNum cos(bfpNum a) { return float_to_bfpNum(std::cos(bfpNum_to_float(a))); }

    bool compare(bfpNum a, bfpNum b)
    {
        return bfpNum_to_float(a) > bfpNum_to_float(b);
        // /* align exponents */
        // if (a.exp >= b.exp)
        // {
        //     b.mant >>= (a.exp - b.exp);
        // }
        // else
        // {
        //     a.mant >>= (b.exp - a.exp);
        // }

        // /* a > b */
        // if (a.sign ^ b.sign)
        // { // 둘 중 하나만 음수인 경우
        //     if (a.sign)
        //         return false; // a가 음수인 경우
        //     else
        //         return true; // b가 음수인 경우
        // }
        // else if (a.sign && b.sign)
        // { // 둘 다 음수인 경우
        //     return a.mant < b.mant;
        // }
        // else
        // { // 둘 다 양수인 경우
        //     return a.mant > b.mant;
        // }
    }

    bool isequal(bfpNum a, bfpNum b)
    {
        return bfpNum_to_float(a) == bfpNum_to_float(b);
        // return (a.sign == b.sign) && (a.exp == b.exp) && (a.mant == b.mant);
    }

    inline bfpNum operator+(bfpNum a, bfpNum b) { return add(a, b); }

    inline bfpNum operator-(bfpNum a, bfpNum b) { return sub(a, b); }

    inline bfpNum operator-(bfpNum a)
    {
        a.sign ^= 1;
        return a;
    }

    inline bfpNum operator*(bfpNum a, bfpNum b) { return mult(a, b); }

    inline bfpNum operator/(bfpNum a, bfpNum b) { return div(a, b); }

    inline bfpNum min(bfpNum a, bfpNum b)
    {
        float a_f = bfpNum_to_float(a);
        float b_f = bfpNum_to_float(b);
        float min = std::fmin(a_f, b_f);
        return float_to_bfpNum(min);
        // /* align exponents */
        // if (a.exp >= b.exp)
        // {
        //     b.mant >>= (a.exp - b.exp);
        // }
        // else
        // {
        //     a.mant >>= (b.exp - a.exp);
        // }

        // /* a > b */
        // if (a.sign ^ b.sign)
        // { // 둘 중 하나만 음수인 경우
        //     if (a.sign)
        //         return a; // a가 음수인 경우
        //     else
        //         return b; // b가 음수인 경우
        // }
        // else if (a.sign && b.sign)
        // { // 둘 다 음수인 경우
        //     if (a.mant < b.mant)
        //         return b;
        //     else
        //         return a;
        // }
        // else
        { // 둘 다 양수인 경우
            if (a.mant > b.mant)
                return b;
            else
                return a;
        }
    }

    inline bfpNum max(bfpNum a, bfpNum b)
    {
        float a_f = bfpNum_to_float(a);
        float b_f = bfpNum_to_float(b);
        float max = std::fmax(a_f, b_f);
        return float_to_bfpNum(max);
        // /* align exponents */
        // if (a.exp >= b.exp)
        // {
        //     b.mant >>= (a.exp - b.exp);
        // }
        // else
        // {
        //     a.mant >>= (b.exp - a.exp);
        // }

        // /* a > b */
        // if (a.sign ^ b.sign)
        // { // 둘 중 하나만 음수인 경우
        //     if (a.sign)
        //         return b; // a가 음수인 경우
        //     else
        //         return a; // b가 음수인 경우
        // }
        // else if (a.sign && b.sign)
        // { // 둘 다 음수인 경우
        //     if (a.mant < b.mant)
        //         return a;
        //     else
        //         return b;
        // }
        // else
        { // 둘 다 양수인 경우
            if (a.mant > b.mant)
                return a;
            else
                return b;
        }
    }

    inline bool operator>(bfpNum a, bfpNum b) { return compare(a, b); }

    inline bool operator>=(bfpNum a, bfpNum b)
    {
        return compare(a, b) || isequal(a, b);
    }

    inline bool operator<(bfpNum a, bfpNum b) { return compare(b, a); }

    inline bool operator<=(bfpNum a, bfpNum b)
    {
        return compare(b, a) || isequal(a, b);
    }

    inline bool operator==(bfpNum a, bfpNum b) { return isequal(a, b); }

    inline bool operator!=(bfpNum a, bfpNum b)
    {
        return (a.sign != b.sign) || (a.exp != b.exp) || (a.mant != b.mant);
    }

    /* arithmetic operations for entire block */
    bfpNum add_bfpBlock_b(bfpBlock block)
    {
        bfpNum res = {0, block.common_exp, 0};
        long long int res_mant_temp = 0;

        // 1. converision to 2's complment for negative mantissas
        for (int i = 0; i < block.sign.size(); i++)
        {
            if (block.sign[i])
            {
                block.M[i] = ~block.M[i] + 1;
            }
        }

        // 2. add mantissas
        for (int i = 0; i < block.M.size(); i++)
        {
            res_mant_temp += (long long int)block.M[i];
            // printf("add ");
            // printBit_sint(block.M[i], false);
            // printf(" => ");
            // printBit_sint(res.mant, true);

            /* option 1 */
            // if((res.mant & 0xff000000) > 0){ //when carry is 1: 11.01 x 2^1
            // = 1.101 x 2^2
            //     res.mant >>= 1;
            //     res.exp += 1;
            //     printf("\tCARRY! : ");
            //     printBit_mant(res.mant, true);
            // }
        }

        // 3. convert to  signed magnitude if negative
        if ((res_mant_temp & 0x8000000000000000) ==
            0x8000000000000000)
        { // if MSB is 1(sign = 1)
            res_mant_temp = ~res_mant_temp + 1;
            res.sign = 1;
        }

        // // if addition result is 0
        // if(res_mant_temp == 0){

        // }

        // 4. normalization
        /* option 2 */
        while ((res_mant_temp >> (BFP_MANT_BITSIZE + 1)) >
               0)
        { // when carry is 1: 11.01 x 2^1 = 1.101 x 2^2
            res_mant_temp >>= 1;
            res.exp += 1;
            // printf("\tCARRY! : ");
            // printBit_bfpNum_exp(res.exp, false);
            // printBit_ulong(res_mant_temp, true);
        }
        res.mant = (int)res_mant_temp;
        int implicit_1 = (int)(std::pow(2, BFP_MANT_BITSIZE));
        while ((res.mant & implicit_1) != implicit_1)
        { // 0.101 x 2^2 = 1.01 x 2^1
            res.mant <<= 1;
            res.exp -= 1;
        }

        return res;
    }

    float add_bfpBlock(bfpBlock block)
    {
        return bfpNum_to_float(add_bfpBlock_b(block));
    }

    bfpNum mult_bfpBlock_b(bfpBlock block)
    {
        bfpNum res = {0,
                      (unsigned int)(block.common_exp * block.M.size() -
                                     ((int)std::pow(2, BFP_EXP_BITSIZE - 1) - 1) *
                                         (block.M.size() - 1)),
                      0};

        // 1. assign sign
        for (int i = 0; i < block.sign.size(); i++)
        {
            res.sign ^= block.sign[i];
        }

        unsigned long long res_mant_temp = (unsigned long long)block.M[0];

        for (int i = 1; i < block.M.size(); i++)
        {
            // 2. multiply mantissas
            res_mant_temp *= (unsigned long long)block.M[i];

            // 3 rounding
            int temp = (int)std::pow(2, BFP_MANT_BITSIZE);
            unsigned short last_bit =
                (unsigned short)((res_mant_temp & temp) >> BFP_MANT_BITSIZE);
            temp >>= 1;
            unsigned short ground_bit =
                (unsigned short)((res_mant_temp & temp) >> BFP_MANT_BITSIZE - 1);
            temp >>= 1;
            unsigned short round_bit =
                (unsigned short)((res_mant_temp & temp) >> BFP_MANT_BITSIZE - 2);
            temp -= 1;
            unsigned short sticky_bits = (unsigned short)(res_mant_temp & temp);

            int lsb = (int)std::pow(2, BFP_MANT_BITSIZE);
            if (ground_bit == 1)
            {
                if (round_bit == 0 && sticky_bits == 0)
                { // round to even
                    if (last_bit == 1)
                    {
                        res_mant_temp += lsb;
                    }
                }
                else
                {
                    res_mant_temp += lsb; // round up
                }
            }

            res_mant_temp >>= BFP_MANT_BITSIZE;
        }

        // 3. normalization
        while ((res_mant_temp >> (BFP_MANT_BITSIZE + 1)) >
               0)
        { // when carry is 1: 11.01 x 2^1 = 1.101 x 2^2
            res_mant_temp >>= 1;
            res.exp += 1;
            // printf("\tCARRY! : ");
            // printBit_bfpNum_exp(res.exp, false);
            // printBit_ulong(res_mant_temp, true);
        }
        res.mant = (int)res_mant_temp;
        int implicit_1 = (int)(std::pow(2, BFP_MANT_BITSIZE));
        while ((res.mant & implicit_1) != implicit_1)
        { // 0.101 x 2^2 = 1.01 x 2^1
            res.mant <<= 1;
            res.exp -= 1;
        }

        return res;
    }

    float mult_bfpBlock(bfpBlock block)
    {
        return bfpNum_to_float(mult_bfpBlock_b(block));
    }
} // namespace bfp
