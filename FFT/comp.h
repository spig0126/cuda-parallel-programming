/* basic complex number operations */
#include <math.h>

typedef struct Comp{
    double r, i;
} Comp;

Comp cal_euler(double x){
    Comp res;
    res.r = cos(x);
    res.i = sin(x);
    return res;
}

Comp comp_mult(Comp a, Comp b){
    Comp res;
    res.r = a.r * b.r - a.i * b.i;
    res.i = a.r * b.i + a.i * b.r;
    return res;
}

Comp comp_add(Comp a, Comp b){
    Comp res;   
    res.r = a.r + b.r;
    res.i = a.i + b.i;
    return res;
}

Comp comp_sub(Comp a, Comp b){
    Comp res;   
    res.r = a.r - b.r;
    res.i = a.i - b.i;
    return res;
}

Comp comp_cpy(Comp a, Comp b){
    b.r = a.r;
    b.i = a.i;
    return a;
}