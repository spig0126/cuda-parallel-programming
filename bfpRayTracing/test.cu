#include <stdio.h>
#include "bfp.cuh"

using namespace std;

// vector<float> f = {0.725133 , 0.729737, 0.671793};
// vector<float> f = {0.725133 , 0.729737};
vector<color> c = {{2.25010704994, 2.250, 2.250107049}, {2.213, 2.2136495, 1.2434}, {6.2341, 2.2131459713, 3.840329}};
// vector<color> c = {{2, 2, 2}, {-2, -2, 2}, {3, 3, 3}};
vector<float> f = {2.25010704994, -3.11704993248};

void compareTwoNum_bfp_float(vector<float> f);

void test_add(vector<float> f);
void test_sub(vector<float> f);
void test_mult(vector<float> f);
void test_div(vector<float> f);

void test_add_float_block(vector<float> f);
void test_mult_float_block(vector<float> f);
void test_add_color_block(vector<color> c);
void test_mult_color_block(vector<color> c);

int main(void){
    test_div(f);

    return 0;
}

void compareTwoNum_bfp_float(vector<float> f, bfpNum a, bfpNum b, bfpNum bfp_res, float float_res){
    printf("=========BFP=============\n");
    printf("%f   =>   ", bfpNum_to_float(a));
    printBit_bfpNum(a, true);
    printf("%f   =>   ", bfpNum_to_float(b));
    printBit_bfpNum(b, true);
    printf("=========FLOAT=============\n");
    printBit_float(f[0]);
    printBit_float(f[1]);

    printf("\nBFP result value:\n");
    printf("%f   =>   ", bfpNum_to_float(bfp_res));
    printBit_bfpNum(bfp_res, true);
    printf("actual value:\n");
    printBit_float(float_res);
}

//--------------------two number arithmetic-----------------------------------
void test_add(vector<float> f){
    bfpNum a = float_to_bfpNum(f[0]);
    bfpNum b = float_to_bfpNum(f[1]);

    bfpNum bfp_res = add(a, b);
    float float_res = f[0] + f[1];

    compareTwoNum_bfp_float(f, a, b, bfp_res, float_res);
}

void test_sub(vector<float> f){
    bfpNum a = float_to_bfpNum(f[0]);
    bfpNum b = float_to_bfpNum(f[1]);

    bfpNum bfp_res = sub(a, b);
    float float_res = f[0] - f[1];

    compareTwoNum_bfp_float(f, a, b, bfp_res, float_res);
}

void test_mult(vector<float> f){
    bfpNum a = float_to_bfpNum(f[0]);
    bfpNum b = float_to_bfpNum(f[1]);

    bfpNum bfp_res = mult(a, b);
    float float_res = f[0] * f[1];

    compareTwoNum_bfp_float(f, a, b, bfp_res, float_res);
}

void test_div(vector<float> f){
    bfpNum a = float_to_bfpNum(f[0]);
    bfpNum b = float_to_bfpNum(f[1]);

    bfpNum bfp_res = div(a, b);
    float float_res = f[0] / f[1];

    compareTwoNum_bfp_float(f, a, b, bfp_res, float_res);
}

//--------------------BFP block arithmetic-----------------------------------
void test_add_float_block(vector<float> f){
    bfpBlock block = createBfpBlock(f);
    print_float_block_formatting(f, block);

    float bfp_res = add_bfpBlock(block);
    float float_res = 0;
    for(int i=0; i<f.size(); i++){
        float_res += f[i];
    }

    printf("\nBFP result value:\n");
    printBit_float(bfp_res);
    printf("\nactual value:\n");
    printBit_float(float_res);
}

void test_mult_float_block(vector<float> f){
    bfpBlock block = createBfpBlock(f);
    print_float_block_formatting(f, block);

    float bfp_res = mult_bfpBlock(block);
    float float_res = 1;
    for(int i=0; i<f.size(); i++){
        float_res *= f[i];
    }

    printf("\nBFP result value:\n");
    printBit_float(bfp_res);
    printf("\nactual value:\n");
    printBit_float(float_res);
}

//--------------------BFP color block arithmetic-----------------------------------
void test_add_color_block(vector<color> colors){
    vector<float> f;
    for(const auto& color : colors){
        f.push_back(color[0]);
        f.push_back(color[1]);
        f.push_back(color[2]);
    }
    bfpBlock block = createBfpBlock(f);
    print_float_block_formatting(f, block);

    color bfp_res = add_color_bfpBlock(block);
    float r=0, g=0, b=0;
    for(int i=0; i<c.size(); i++){
        r += c[i][0];
        g += c[i][1];
        b += c[i][2];
    }

    printf("\nBFP result value:\n");
    print_color(bfp_res);
    printf("\nactual value:\n");
    printBit_float(r);
    printBit_float(g);
    printBit_float(b);
    return;
}

void test_mult_color_block(vector<color> colors){
    vector<float> f;
    for(const auto& color : colors){
        f.push_back(color[0]);
        f.push_back(color[1]);
        f.push_back(color[2]);
    }
    bfpBlock block = createBfpBlock(f);
    print_float_block_formatting(f, block);

    color bfp_res = add_color_bfpBlock(block);
    float r=1, g=1, b=1;
    for(int i=0; i<c.size(); i++){
        r *= c[i][0];
        g *= c[i][1];
        b *= c[i][2];
    }

    printf("\nBFP result value:\n");
    print_color(bfp_res);
    printf("\nactual value:\n");
    printBit_float(r);
    printBit_float(g);
    printBit_float(b);
    return;
}