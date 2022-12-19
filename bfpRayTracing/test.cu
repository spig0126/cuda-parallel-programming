#include <stdio.h>
#include "bfp.cuh"

vector<float> f = {0.725133 , 0.729737, 0.671793};
// vector<color> c = {{2.25010704994, 2.250, 2.250107049}, {2.213, 2.2136495, 1.2434}, {6.2341, 2.2131459713, 3.840329}};
// vector<color> c = {{2, 2, 2}, {-2, -2, 2}, {3, 3, 3}};


int main(void){
    test_add_float_block(f);
    // test_mult_float_block(f);

    // float f1 =0.725133;
    // bfpNumFloat bfp_f1 = float_to_bfpNumFloat(f1);
    // printBit_sint(bfp_f1.mant, true);

    // float f2 =0.729737 + f1;
    // bfpNumFloat bfp_f2 = float_to_bfpNumFloat(f2);
    // printBit_sint(bfp_f2.mant, true);
    
    // float f3 =0.671793 + f2;
    // bfpNumFloat bfp_f3 = float_to_bfpNumFloat(f3);
    // printBit_sint(bfp_f3.mant, true);
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