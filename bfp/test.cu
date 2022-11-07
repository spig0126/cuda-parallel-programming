#include <stdio.h>
#include "bfp.cuh"

float a[N] = {2, 3};


void two_num_add(bfpBlock block){
    bfpNumFloat res = add_f(block, block.M[0], block.M[1]);
    printBit_bfpNumFloat(res, true);
    print_bfpNumFloat(res);
    float res_float = bfpNumFloat_to_float(res);
    printf("%f\n", res_float);
}

void block_add(bfpBlock block){
    bfpNumFloat res = add_block_f(block);
    float res_float = bfpNumFloat_to_float(res);
    printf("\n\n=====================RESULT============================\nbit representation: ");
    printBit_bfpNumFloat(res, true);
    printf("\n");
    print_bfpNumFloat(res);
    printf("\nvalue = %f\n", res_float);
}


void two_num_mult(bfpBlock block){
    bfpNumFloat res = mult_f(block, block.M[0], block.M[1]);
    float res_float = bfpNumFloat_to_float(res);
    printf("\n\n=====================RESULT============================\n");
    printf("value = %.9g\n", res_float);
    printf("bit representation: ");
    printBit_bfpNumFloat(res, true);
    printf("\n");
    printf("<FP multiplication>\nvalue = %.9g\nbit representation: ", a[0] * a[1]);
    printBit_float(a[0] * a[1]);
    printf("\n");
    print_bfpNumFloat(res);

}

void two_num_div(bfpBlock block){
    bfpNumFloat res = div_f(block, block.M[0], block.M[1]);
    float res_float = bfpNumFloat_to_float(res);
    printf("\n\n=====================RESULT============================\n");
    printf("value = %.9g\n", res_float);
    printf("bit representation: ");
    printBit_bfpNumFloat(res, true);
    printf("\n");
    printf("<FP division>\nvalue = %.9g\nbit representation: ", a[0] / a[1]);
    printBit_float(a[0] / a[1]);
    printf("\n");
    print_bfpNumFloat(res);
}

int main(void){
    /* block formatting */
    printf("---------Before block formatiing---------\n");
    for(int i=0; i<N; i++){
        printBit_float(a[i]);
    }
    bfpBlock block = createBfpBlock(a);
    print_bfpBlock(block);
    printf("===============================================\n");

    two_num_mult(block);
    return 0;
}