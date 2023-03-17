#ifndef SPHERE_BFP_H
#define SPHERE_BFP_H

#include "hittable_bfp.h"

class sphere : public hittable
{
public:
    sphere() {}
    sphere(int idx, point3 cen, bfpNum r, shared_ptr<material> m) : sphereIdx(idx), center(cen), radius(r), mat_ptr(m){};

    virtual bool hit(const ray &r, bfpNum t_min, bfpNum t_max, hit_record &rec) const override;
    virtual bool bounding_box(bfpNum time0, bfpNum time1, aabb &output_box) const override;

    point3_float center_f() const { return vec3_bfpNum_to_float(center); }
    float radius_f() const { return bfpNum_to_float(radius); }

public:
    point3 center;
    bfpNum radius;
    shared_ptr<material> mat_ptr;
    int sphereIdx;
};

bool sphere::hit(const ray &r, bfpNum t_min, bfpNum t_max, hit_record &rec) const
{
    /* before: bfpNUm */
    vec3 oc = r.origin() - center;
    bfpNum a = r.direction().length_squared();
    bfpNum half_b = dot(oc, r.direction());
    bfpNum c = oc.length_squared() - radius * radius;
    bfpNum discriminant = half_b * half_b - a * c;
    bfpNum sqrtd = sqrt(discriminant);
    bfpNum root = (-half_b - sqrtd) / a;

    point3 pm = center - vec3(radius, radius, radius);
    point3 pM = center + vec3(radius, radius, radius);

    if (discriminant < b_0) // if the ray doesn't hit the sphere
        return false;
    if (root < t_min || t_max < root) // If the ray hits the sphere,
    {
        root = (-half_b + sqrtd) / a; // Find the nearest root that lies in the acceptable range.
        if (root < t_min || t_max < root)
            return false;
    }

    rec.t = root;
    rec.p = r.at(rec.t);
    vec3 outward_normal = (rec.p - center) / radius;
    rec.set_face_normal(r, outward_normal);
    rec.mat_ptr = mat_ptr;

    return true;

    // /* after: float */
    // float t_min_f = bfpNum_to_float(t_min);
    // float t_max_f = bfpNum_to_float(t_max);

    // vec3_float oc = r.origin_f() - center_f();
    // float a = r.direction_f().length_squared();
    // float half_b = dot(oc, r.direction_f());
    // float c = oc.length_squared() - radius_f() * radius_f();
    // float discriminant = half_b * half_b - a * c;
    // float sqrt_d = std::sqrt(discriminant);
    // float root = (-half_b - sqrt_d) / a;

    // point3_float pm = center_f() - vec3_float(radius_f(), radius_f(), radius_f());
    // point3_float pM = center_f() + vec3_float(radius_f(), radius_f(), radius_f());

    // if (discriminant < 0.0)
    //     return false;
    // if (root < t_min_f || t_max_f < root)
    // {
    //     root = (-half_b + sqrt_d) / a;
    //     if (root < t_min_f || t_max_f < root)
    //         return false;
    // }

    // rec.t = float_to_bfpNum(root);
    // rec.p = r.at(rec.t);

    // vec3_float p_f = vec3_bfpNum_to_float(rec.p);
    // vec3_float outward_normal = (p_f - center_f()) / radius_f();
    // rec.set_face_normal(r, vec3_float_to_bfpNum(outward_normal));
    // rec.mat_ptr = mat_ptr;

    // return true;
}

bool sphere::bounding_box(bfpNum time0, bfpNum time1, aabb &output_box) const
{

    output_box = aabb(
        center - vec3(radius, radius, radius),
        center + vec3(radius, radius, radius));

    return true;
}

#endif
