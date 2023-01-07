#ifndef MOVING_SPHERE_BFP_H
#define MOVING_SPHERE_BFP_H

#include "hittable_bfp.h"

class moving_sphere : public hittable
{
public:
    moving_sphere() {}
    moving_sphere(
        int idx, point3 cen0, point3 cen1, bfpNum _time0, bfpNum _time1, bfpNum r, shared_ptr<material> m)
        : sphereIdx(idx), center0(cen0), center1(cen1), time0(_time0), time1(_time1), radius(r), mat_ptr(m){};

    virtual bool hit(const ray &r, bfpNum t_min, bfpNum t_max, hit_record &rec) const override;
    virtual bool bounding_box(bfpNum _time0, bfpNum _time1, aabb &output_box) const override;
    point3 center(bfpNum time) const;

public:
    point3 center0, center1;
    bfpNum time0, time1;
    bfpNum radius;
    shared_ptr<material> mat_ptr;
    int sphereIdx;
};

point3 moving_sphere::center(bfpNum time) const
{
    return center0 + ((time - time0) / (time1 - time0)) * (center1 - center0);
}

bool moving_sphere::hit(const ray &r, bfpNum t_min, bfpNum t_max, hit_record &rec) const
{
    vec3 oc = r.origin() - center(r.time());
    auto a = r.direction().length_squared();
    auto half_b = dot(oc, r.direction());
    auto c = oc.length_squared() - radius * radius;

    auto discriminant = half_b * half_b - a * c;
    if (discriminant < b_0)
        return false;
    auto sqrtd = sqrt(discriminant);

    // Find the nearest root that lies in the acceptable range.
    auto root = (-half_b - sqrtd) / a;
    if (root < t_min || t_max < root)
    {
        root = (-half_b + sqrtd) / a;
        if (root < t_min || t_max < root)
            return false;
    }

    rec.t = root;
    rec.p = r.at(rec.t);
    auto outward_normal = (rec.p - center(r.time())) / radius;
    rec.set_face_normal(r, outward_normal);
    rec.mat_ptr = mat_ptr;

    return true;
}

bool moving_sphere::bounding_box(bfpNum _time0, bfpNum _time1, aabb &output_box) const
{
    aabb box0(
        center(_time0) - vec3(radius, radius, radius),
        center(_time0) + vec3(radius, radius, radius));
    aabb box1(
        center(_time1) - vec3(radius, radius, radius),
        center(_time1) + vec3(radius, radius, radius));
    output_box = surrounding_box(box0, box1);
    return true;
}

#endif
