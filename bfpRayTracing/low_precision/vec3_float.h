/*
 * ===================================================
 *
 *       Filename:  vec3.h
 *    Description:  Ray Tracing: The Next Week (RTTNW): ~BVH
 *        Created:  2022/07/13`
 *
 * ===================================================
 */

// Preprocessors

#ifndef VEC3_FLOAT_H
#define VEC3_FLOAT_H
#include <cmath>
#include <iostream>
#include "utility_bfp.h"

// Usings

// Classes

class vec3_float
{
public:
    vec3_float() : e{0, 0, 0} {}
    vec3_float(float e0, float e1, float e2) : e{e0, e1, e2} {}

    float x() const { return e[0]; }
    float y() const { return e[1]; }
    float z() const { return e[2]; }

    vec3_float operator-() const { return vec3_float(-e[0], -e[1], -e[2]); }
    float operator[](int i) const { return e[i]; }
    float &operator[](int i) { return e[i]; }

    vec3_float &operator+=(const vec3_float &v)
    {
        e[0] += v.e[0];
        e[1] += v.e[1];
        e[2] += v.e[2];

        return *this;
    }

    vec3_float &operator*=(const float t)
    {
        e[0] *= t;
        e[1] *= t;
        e[2] *= t;

        return *this;
    }

    vec3_float &operator/=(const float t)
    {
        return *this *= 1 / t;
    }

    float length() const
    {
        return std::sqrt(length_squared());
    }

    float length_squared() const
    {
        return e[0] * e[0] + e[1] * e[1] + e[2] * e[2];
    }

    bool near_zero(bool t1, bool t2) const
    {
        // Return true if the vector is close to zero in all dimensions.
        const auto s = 1e-8;
        return (fabs(e[0]) < s) && (fabs(e[1]) < s) && (fabs(e[2]) < s);
    }

    static vec3_float random()
    {
        return vec3_float(random_float(), random_float(), random_float());
    }

    static vec3_float random(float min, float max)
    {
        return vec3_float(random_float(min, max), random_float(min, max), random_float(min, max));
    }

public:
    float e[3];
};

// Utility Functions

inline std::ostream &operator<<(std::ostream &out, const vec3_float &v)
{
    return out << v.e[0] << ' ' << v.e[1] << ' ' << v.e[2];
}

inline vec3_float operator+(const vec3_float &u, const vec3_float &v)
{
    return vec3_float(u.e[0] + v.e[0], u.e[1] + v.e[1], u.e[2] + v.e[2]);
}

inline vec3_float operator-(const vec3_float &u, const vec3_float &v)
{
    return vec3_float(u.e[0] - v.e[0], u.e[1] - v.e[1], u.e[2] - v.e[2]);
}

inline vec3_float operator*(const vec3_float &u, const vec3_float &v)
{
    return vec3_float(u.e[0] * v.e[0], u.e[1] * v.e[1], u.e[2] * v.e[2]);
}

inline vec3_float operator*(float t, const vec3_float &v)
{
    return vec3_float(t * v.e[0], t * v.e[1], t * v.e[2]);
}

inline vec3_float operator*(const vec3_float &v, float t)
{
    return t * v;
}

inline vec3_float operator/(vec3_float v, float t)
{
    return (1 / t) * v;
}

inline float dot(const vec3_float &u, const vec3_float &v)
{
    return u.e[0] * v.e[0] + u.e[1] * v.e[1] + u.e[2] * v.e[2];
}

inline vec3_float cross(const vec3_float &u, const vec3_float &v)
{
    return vec3_float(u.e[1] * v.e[2] - u.e[2] * v.e[1],
                      u.e[2] * v.e[0] - u.e[0] * v.e[2],
                      u.e[0] * v.e[1] - u.e[1] * v.e[0]);
}

inline vec3_float unit_vector(vec3_float v)
{
    return v / v.length();
}

inline vec3_float random_in_unit_sphere(bool t1, bool t2)
{
    while (true)
    {
        vec3_float p = vec3_float::random(-1, 1);
        if (p.length_squared() >= 1)
            continue;

        return p;
    }
}

inline vec3_float random_unit_vector(bool t1, bool t2)
{
    return unit_vector(random_in_unit_sphere(true, true));
}

inline vec3_float random_in_hemisphere(const vec3_float &normal)
{
    vec3_float in_unit_sphere = random_in_unit_sphere(true, true);

    if (dot(in_unit_sphere, normal) > 0.0) // In the same hemisphere as the normal
        return in_unit_sphere;
    else
        return -in_unit_sphere;
}

inline vec3_float reflect(const vec3_float &v, const vec3_float &n)
{
    return v - 2 * dot(v, n) * n;
}

inline vec3_float refract(const vec3_float &uv, const vec3_float &n, float etai_over_etat)
{
    auto cos_theta = fmin(dot(-uv, n), 1.0);
    vec3_float r_out_perp = etai_over_etat * (uv + cos_theta * n);
    vec3_float r_out_parallel = -sqrt(fabs(1.0 - r_out_perp.length_squared())) * n;
    return r_out_perp + r_out_parallel;
}

vec3_float random_in_unit_disk()
{
    while (true)
    {
        auto p = vec3_float(random_float(-1, 1), random_float(-1, 1), 0);
        if (p.length_squared() >= 1)
            continue;
        return p;
    }
}

// Type aliases for vec3_float

using point3_float = vec3_float; // 3D point
using color_float = vec3_float;  // RGB color

#endif
