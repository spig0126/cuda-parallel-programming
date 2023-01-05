#ifndef MATERIAL_BFP_H
#define MATERIAL_BFP_H

#include "hittable_bfp.h"

using namespace bfp;

class material
{
public:
    virtual bool scatter(const ray &r_in, const hit_record &rec, color &attenuation, ray &scattered) const = 0;
};

class lambertian : public material
{
public:
    lambertian(const color &a) : albedo(a) {}

    virtual bool scatter(const ray &r_in, const hit_record &rec, color &attenuation, ray &scattered) const override
    {
        vec3 scatter_direction = rec.normal + random_unit_vector(true);

        // Catch degenerate scatter direction
        if (scatter_direction.near_zero())
            scatter_direction = rec.normal;

        scattered = ray(rec.p, scatter_direction, r_in.time());
        attenuation = albedo;
        return true;
    }

public:
    color albedo;
};

class metal : public material
{
public:
    metal(const color &a, bfpNum f) : albedo(a), fuzz(f < b_1 ? f : b_1) {}

    virtual bool scatter(const ray &r_in, const hit_record &rec, color &attenuation, ray &scattered) const override
    {
        vec3 reflected = reflect(unit_vector(r_in.direction()), rec.normal);
        scattered = ray(rec.p, reflected + fuzz * random_in_unit_sphere(true), r_in.time());
        attenuation = albedo;
        return (dot(scattered.direction(), rec.normal) > b_0);
    }

public:
    color albedo;
    bfpNum fuzz;
};

class dielectric : public material
{
public:
    dielectric(bfpNum index_of_refraction) : ir(index_of_refraction) {}

    virtual bool scatter(
        const ray &r_in, const hit_record &rec, color &attenuation, ray &scattered) const override
    {
        attenuation = color(b_1, b_1, b_1);
        bfpNum refraction_ratio = rec.front_face ? (b_1 / ir) : ir;

        vec3 unit_direction = unit_vector(r_in.direction());

        bfpNum cos_theta = min(dot(-unit_direction, rec.normal), b_1);
        bfpNum sin_theta = sqrt(b_1 - cos_theta * cos_theta);

        bool cannot_refract = refraction_ratio * sin_theta > b_1;
        vec3 direction;

        if (cannot_refract || reflectance(cos_theta, refraction_ratio) > random_num())
            direction = reflect(unit_direction, rec.normal);
        else
            direction = refract(unit_direction, rec.normal, refraction_ratio);

        scattered = ray(rec.p, direction, r_in.time());
        return true;
    }

public:
    bfpNum ir;

private:
    static bfpNum reflectance(bfpNum cosine, bfpNum ref_idx)
    {
        // Use Schlick's approximation for reflectance.
        bfpNum r0 = (b_1 - ref_idx) / (b_1 + ref_idx);
        r0 = r0 * r0;
        return r0 + (b_1 - r0) * pow((b_1 - cosine), 5);
    }
};

#endif