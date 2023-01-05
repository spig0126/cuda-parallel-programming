#ifndef HITTABLE_BFP_H
#define HITTABLE_BFP_H

#include <memory>

#include "aabb_bfp.h"
#include "utility_bfp.h"

using namespace bfp;
using std::make_shared;
using std::shared_ptr;

class material;

struct hit_record
{
    point3 p;
    vec3 normal;
    shared_ptr<material> mat_ptr;
    bfpNum t;
    bool front_face;

    inline void set_face_normal(const ray &r, const vec3 &outward_normal)
    {
        front_face = dot(r.direction(), outward_normal) < b_0;
        normal = front_face ? outward_normal : -outward_normal;
    }
};

class hittable
{
public:
    virtual bool hit(const ray &r, bfpNum t_min, bfpNum t_max, hit_record &rec) const = 0;
    virtual bool bounding_box(bfpNum time0, bfpNum time1, aabb &output_box) const = 0;
};

#endif