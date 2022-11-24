#include <stdio.h>
#include "bfp.cuh"

vector<float> a = {2.2, 2.2};
vector<color> c = {{2.25010704994, 2.25010704994, 2.25010704994}, {2.21364951134, 2.21364951134, 2.21364951134}, {2.2131459713, 2.2131459713, 2.2131459713}};
// vector<color> c = {{2, 2, 2}, {-2, -2, 2}, {3, 3, 3}};



int main(void){
    /* block formatting */
    printf("---------Before block formatiing---------\n");
    for(int i=0; i<c.size(); i++){
        printBit_float(c[i][0]);
        printBit_float(c[i][1]);
        printBit_float(c[i][2]);
    }
    bfpBlock block = createColorBfpBlock(c);
    print_bfpBlock(block);
    printf("===============================================\n");

    color res = mult_color_bfpBlock(block);
    print_color(res);
    printf("\nactual values:\n");
    float r=1, g=1, b=1;
    for(int i=0; i<c.size(); i++){
        r *= c[i][0];
        g *= c[i][1];
        b *= c[i][2];
    }
    printBit_float(r);
    printBit_float(g);
    printBit_float(b);
    return 0;
}


// void two_num_add(bfpBlock block){
//     bfpNumFloat res = add_f(block, block.M[0], block.M[1]);
//     printBit_bfpNumFloat(res, true);
//     print_bfpNumFloat(res);
//     float res_float = bfpNumFloat_to_float(res);
//     printf("%f\n", res_float);
// }

// void block_add(bfpBlock block){
//     bfpNumFloat res = add_block_f(block);
//     float res_float = bfpNumFloat_to_float(res);
//     printf("\n\n=====================RESULT============================\nbit representation: ");
//     printBit_bfpNumFloat(res, true);
//     printf("\n");
//     print_bfpNumFloat(res);
//     printf("\nvalue = %f\n", res_float);
// }


// void two_num_mult(bfpBlock block){
//     bfpNumFloat res = mult_f(block, block.M[0], block.M[1]);
//     float res_float = bfpNumFloat_to_float(res);
//     printf("\n\n=====================RESULT============================\n");
//     printf("value = %.9g\n", res_float);
//     printf("bit representation: ");
//     printBit_bfpNumFloat(res, true);
//     printf("\n");
//     printf("<FP multiplication>\nvalue = %.9g\nbit representation: ", a[0] * a[1]);
//     printBit_float(a[0] * a[1]);
//     printf("\n");
//     print_bfpNumFloat(res);

// }

// void two_num_div(bfpBlock block){
//     bfpNumFloat res = div_f(block, block.M[0], block.M[1]);
//     float res_float = bfpNumFloat_to_float(res);
//     printf("\n\n=====================RESULT============================\n");
//     printf("value = %.9g\n", res_float);
//     printf("bit representation: ");
//     printBit_bfpNumFloat(res, true);
//     printf("\n");
//     printf("<FP division>\nvalue = %.9g\nbit representation: ", a[0] / a[1]);
//     printBit_float(a[0] / a[1]);
//     printf("\n");
//     print_bfpNumFloat(res);
// }