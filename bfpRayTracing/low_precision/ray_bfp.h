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

    point3 at(bfpNum t) const
    {
        return orig + t * dir;
    }

public:
    point3 orig;
    vec3 dir;
    bfpNum tm;
};

#endif