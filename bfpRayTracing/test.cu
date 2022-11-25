#include <stdio.h>
#include "bfp.cuh"

vector<float> f = {2.25010704994, 2.250, 2.250107049, 2.213, 2.2136495, 1.2434, 6.2341, 2.2131459713, 3.840329};
// vector<color> c = {{2.25010704994, 2.250, 2.250107049}, {2.213, 2.2136495, 1.2434}, {6.2341, 2.2131459713, 3.840329}};
// vector<color> c = {{2, 2, 2}, {-2, -2, 2}, {3, 3, 3}};
// vector<float> f = {2.25010704994, -3.11704993248};

void test_add_float_block(vector<float> f);
void test_mult_float_block(vector<float> f);
void test_add_color_block(vector<color> c);
void test_mult_color_block(vector<color> c);

int main(void){
    test_add_float_block(f);
    return 0;
}

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

void test_add_color_block(vector<color> c){
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

void test_mult_color_block(vector<color> c){
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