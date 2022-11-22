#include <stdio.h>
#include <cuComplex.h>
#include <type_traits>
#include <typeinfo>
#include <stdbool.h>

#define isCompatible(x, type) _Generic(x, type: true, default: false)

int main(void){
    typedef float fp;

    printf("%d\n", sizeof(fp));

}