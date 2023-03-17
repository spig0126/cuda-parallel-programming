#pragma once
#ifndef RAY_BFP_H
#define RAY_BFP_H

#include "vec3_bfp.h"

using namespace bfp;

class ray
{
public:
    ray() {}
    ray(const point3 &origin, const point3 &direction, bfpNum time) : orig(origin), dir(direction), tm(time) {}

    point3 origin() const { return orig; }
    vec3 direction() const { return dir; }
    bfpNum time() const { return tm; }

    point3_float origin_f() const { return vec3_bfpNum_to_float(orig); }
    vec3_float direction_f() const { return vec3_bfpNum_to_float(dir); }
    float time_f() const { return bfpNum_to_float(tm); }

    point3 at(bfpNum t) const
    {
        return orig + t * dir;
    }

    point3_float at(float t) const
    {
        return origin_f() + t * direction_f();
    }

public:
    point3 orig;
    vec3 dir;
    bfpNum tm;
};

#endif