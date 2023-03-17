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
        /* before (calculated with bfpNum) */
        vec3 scatter_direction = rec.normal + random_unit_vector(true);
        // Catch degenerate scatter direction
        if (scatter_direction.near_zero())
            scatter_direction = rec.normal;

        scattered = ray(rec.p, scatter_direction, r_in.time());
        attenuation = albedo;
        return true;

        // /* after (calculated with float) */
        // vec3_float scatter_direction_f = vec3_bfpNum_to_float(rec.normal) + random_unit_vector(true, true);
        // vec3 scatter_direction = vec3_float_to_bfpNum(scatter_direction_f);

        // // Catch degenerate scatter direction
        // if (scatter_direction.near_zero())
        //     scatter_direction = rec.normal;

        // scattered = ray(rec.p, scatter_direction, r_in.time());
        // attenuation = albedo;
        // return true;
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
        /* before (calculated with bfpNum) */
        vec3 reflected = reflect(unit_vector(r_in.direction()), rec.normal);
        scattered = ray(rec.p, reflected + fuzz * random_in_unit_sphere(true), r_in.time());
        attenuation = albedo;
        return (dot(scattered.direction(), rec.normal) > b_0);

        // /* after (calculated with float) */
        // vec3_float reflected = reflect(unit_vector(r_in.direction_f()), vec3_bfpNum_to_float(rec.normal));
        // scattered = ray(rec.p, vec3_float_to_bfpNum(reflected + bfpNum_to_float(fuzz) * random_in_unit_sphere(true, true)), r_in.time());
        // attenuation = albedo;
        // return (dot(scattered.direction_f(), vec3_bfpNum_to_float(rec.normal)) > 0.0);
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
        /* before (calculated with bfpNum) */
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

        // /* after */
        // float ir_f = bfpNum_to_float(ir);
        // vec3_float normal_f = vec3_bfpNum_to_float(rec.normal);

        // attenuation = color(b_1, b_1, b_1);
        // float refraction_ratio = rec.front_face ? (1.0 / ir_f) : ir_f;

        // vec3_float unit_direction = unit_vector(r_in.direction_f());
        // float cos_theta = fmin(dot(-unit_direction, normal_f), 1.0);
        // float sin_theta = sqrtf(1.0 - cos_theta * cos_theta);

        // bool cannot_refract = refraction_ratio * sin_theta > 1.0;
        // vec3_float direction;

        // if (cannot_refract || reflectance(cos_theta, refraction_ratio) > random_float())
        //     direction = reflect(unit_direction, normal_f);
        // else
        //     direction = refract(unit_direction, normal_f, refraction_ratio);

        // scattered = ray(rec.p, vec3_float_to_bfpNum(direction), r_in.time());
        // return true;
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

    static float reflectance(float cosine, float ref_idx)
    {
        // Use Schlick's approximation for reflectance.
        auto r0 = (1 - ref_idx) / (1 + ref_idx);
        r0 = r0 * r0;
        return r0 + (1 - r0) * pow((1 - cosine), 5);
    }
};

#endif