#include "bfpStruct.cuh"
#include "bfp.cuh"

class vec3_bfp{
    bfpNum e[3];

    public:
        vec3_bfp() : e { {0, 0, 0}, {0, 0, 0}, {0, 0, 0}} {}
        vec3_bfp(bfpNum e0, bfpNum e1, bfpNum e2) : e{e0, e1, e2} {}

        bfpNum x() const { return e[0]; }
        bfpNum y() const { return e[1]; }
        bfpNum z() const { return e[2]; }
        
}