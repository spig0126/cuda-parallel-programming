#pragma once
#include <string.h>
#include <iostream>

#include "bfp.cuh"

void print_bfpNum(bfpNum b);
void print_bfpBlock(bfpBlock block);
void printBit_bfpNum(bfpNum b, bool nextLine);
void printBit_bfpNum_exp(unsigned int e, bool nextLine);
void printBit_bfpNum_mant(int m, bool nextLine);
void printBit_float_mant(int m, bool newLine);
void printBit_uint(unsigned int num, int len, bool reverse);
void printBit_float(float f);
void printBit_sint(int num, bool newLine);
void printBit_ulong(long long num, bool newLine);

using namespace std;

/* print Structs*/
void print_bfpNum(bfpNum b){
    cout << "---------BFP component -------------" << endl;
    cout << "sign: " << b.sign << endl;
    cout << "exp: " << b.exp << endl;
    cout << "mant: " << b.mant << endl;
}

void print_bfpBlock(bfpBlock block){
    cout << "----------BFP block ------------------\n";
    cout << "N: " << block.M.size() << "\ncommon_exp: ";
    printBit_uint(block.common_exp, FLOAT_EXP_BITSIZE, false);
    cout << endl;
    for(int i=0; i<block.M.size(); i++){
        cout << i << ": ";
        printBit_uint(block.sign[i], 1, false);
        cout << "\t";
        printBit_bfpNum_mant(block.M[i], true);
    }
}

/* print processes */
void print_float_block_formatting(vector<float> f, bfpBlock block){
    cout << "=========================================================" << endl;
    cout << "---------Before block formatiing---------" << endl;

    for(int i=0; i<f.size(); i++){
        printBit_float(f[i]);
    }
    print_bfpBlock(block);
    cout << "=========================================================" << endl;
}

/* print bit representations of bfpStructs */
void printBit_bfpNum(bfpNum b, bool nextLine){
    printBit_uint(b.sign, 1, false);
    cout << " " ;
    printBit_uint(b.exp, BFP_EXP_BITSIZE, false);
    cout << " " ;
    printBit_uint(b.mant, BFP_MANT_BITSIZE, false);

    if(nextLine){
        cout << endl;
    }
}

void printBit_bfpNum_exp(unsigned int e, bool nextLine){
    cout << e << " => ";
    printBit_uint(e, BFP_EXP_BITSIZE, false);
    
    if(nextLine){
        cout << endl;
    }
}

void printBit_bfpNum_mant(int m, bool nextLine){
    int temp = m;
    vector<string> out;
    for(int i=0; i<BFP_MANT_BITSIZE + 1; i++, temp>>=1){
        if(i%4 == 0){
            out.insert(out.begin(), " ");
        }
        if(temp&1){
            out.insert(out.begin(), "1");
        }
        else{
            out.insert(out.begin(), "0");
        }
    }

    for(const auto& bit: out){
        cout << bit;
    }

    if(nextLine){
        cout << "\n";
    }
}

void printBit_float_mant(int m, bool newLine){
    int temp = m;
    char out[MAX_BITSIZE] = "";
    for(int i=0; i<24; i++, temp>>=1){
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
    
    for(int i=29; i>=0; i--){
        cout << out[i];
    }
    printf(" (0x%.8x)", m);

    if(newLine){
        cout << "\n";
    }
}

/* print bit representations of basic data types */
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
            cout << out[i];
        }
        
    }
    else{
        for(int i=len-1; i>=0; i--){
            cout << out[i];
        }
    }
}

void printBit_float(float f){
    cout << f << "   =>   " ;
    //cast to integer for bitwise operations
    unsigned int* temp = reinterpret_cast<unsigned int*>(&f); 
    unsigned int orginal_val = *temp;

    char out[MAX_BITSIZE] = "";
    for(int i=0; i<32; i++, *temp>>=1){
        if(i == 23 || i == 31){
            strcat(out, " ");
        }
        if(*temp&1){
            strcat(out, "1");
        }
        else{
            strcat(out, "0");
        }
    }
    
    for(int i=35; i>=0; i--){
            cout << out[i];
    }
    printf(" (0x%.8x)\n", orginal_val);
}

void printBit_sint(int num, bool newLine){     //4bits grouped together
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
        cout << out[i];
    }   
    printf(" (0x%x)", num);

    if(newLine){
        cout << endl;
    }
}

void printBit_ulong(long long num, bool newLine){   //4bits grouped together
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
        cout << out[i];
    }

    printf("\t(%llx)", num);

    if(newLine){
        cout << endl;
    }
}
