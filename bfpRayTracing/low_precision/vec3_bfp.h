#pragma once
#ifndef VEC3_BFP_H
#define VEC3_BFP_H

#include "utility_bfp.h"
#include <iostream>

using namespace bfp;

class vec3
{
public:
    vec3() : e{{0, 0, 0}, {0, 0, 0}, {0, 0, 0}} {}
    vec3(bfpNum e0, bfpNum e1, bfpNum e2) : e{e0, e1, e2} {}

    bfpNum x() const { return e[0]; }
    bfpNum y() const { return e[1]; }
    bfpNum z() const { return e[2]; }

    vec3 operator-() const { return vec3(abs(e[0]), abs(e[1]), abs(e[2])); }
    bfpNum operator[](int i) const { return e[i]; }
    bfpNum &operator[](int i) { return e[i]; }

    vec3 &operator+=(const vec3 &v)
    {
        e[0] = e[0] + v.e[0];
        e[1] = e[1] + v.e[1];
        e[2] = e[2] + v.e[2];

        return *this;
    }

    vec3 &operator*=(const bfpNum t)
    {
        e[0] = e[0] * t;
        e[1] = e[1] * t;
        e[2] = e[2] * t;

        return *this;
    }

    vec3 &operator/=(bfpNum t)
    {
        bfpNum temp = b_1 / t;
        *this = vec3(e[0] * temp, e[1] * temp, e[2] * temp);
        return *this;
    }

    bfpNum length() const
    {
        return sqrt(length_squared());
    }
    bfpNum length_squared() const
    {
        return e[0] * e[0] + e[1] * e[1] + e[2] * e[2];
    }

    bool near_zero() const
    {
        bfpNum s = int_to_bfpNum(0x322bcc77);
        return (abs(e[0]) < s) && (abs(e[1]) < s) && (abs(e[2]) < s);
    }

    static vec3 random()
    {
        return vec3(random_num(), random_num(), random_num());
    }

    static vec3 random(bfpNum min, bfpNum max)
    {
        return vec3(random_num(min, max), random_num(min, max), random_num(min, max));
    }

public:
    bfpNum e[3];
};

// utility functions
inline std::ostream &operator<<(std::ostream &out, const vec3 &v)
{
    return out << bfpNum_to_float(v.e[0]) << ' ' << bfpNum_to_float(v.e[1]) << ' ' << bfpNum_to_float(v.e[2]);
}

inline vec3 operator+(const vec3 &u, const vec3 &v)
{
    return vec3(u.e[0] + v.e[0], u.e[1] + v.e[1], u.e[2] + v.e[2]);
}

inline vec3 operator-(const vec3 &u, const vec3 &v)
{
    return vec3(u.e[0] - v.e[0], u.e[1] - v.e[1], u.e[2] - v.e[2]);
}

inline vec3 operator*(const vec3 &u, const vec3 &v)
{
    return vec3(u.e[0] * v.e[0], u.e[1] * v.e[1], u.e[2] * v.e[2]);
}

inline vec3 operator*(bfpNum t, const vec3 &v)
{
    return vec3(t * v.e[0], t * v.e[1], t * v.e[2]);
}

inline vec3 operator*(const vec3 &v, bfpNum t)
{
    return t * v;
}

inline vec3 operator/(vec3 v, bfpNum t)
{
    return (b_1 / t) * v;
}

inline bool isequal(vec3 v1, vec3 v2)
{
    return v1[0] == v2[0] && v1[1] == v2[1] && v1[2] == v2[2];
}

inline bfpNum dot(const vec3 &u, const vec3 &v)
{
    return u.e[0] * v.e[0] + u.e[1] * v.e[1] + u.e[2] * v.e[2];
}

inline vec3 cross(const vec3 &u, const vec3 &v)
{
    return vec3(u.e[1] * v.e[2] - u.e[2] * v.e[1],
                u.e[2] * v.e[0] - u.e[0] * v.e[2],
                u.e[0] * v.e[1] - u.e[1] * v.e[0]);
}

inline vec3 unit_vector(vec3 v)
{
    return v / v.length();
}

inline vec3 random_in_unit_sphere(bool t)
{
    while (true)
    {
        vec3 p = vec3::random(b_1_neg, b_1);
        if (p.length_squared() >= b_1)
            continue;

        return p;
    }
}

inline vec3 random_unit_vector(bool t)
{
    return unit_vector(random_in_unit_sphere(true));
}

inline vec3 random_in_hemisphere(const vec3 &normal)
{
    vec3 in_unit_sphere = random_in_unit_sphere(true);

    if (dot(in_unit_sphere, normal) > b_0) // In the same hemisphere as the normal
        return in_unit_sphere;
    else
        return -in_unit_sphere;
}

inline vec3 reflect(const vec3 &v, const vec3 &n)
{
    return v - b_2 * dot(v, n) * n;
}

inline vec3 refract(const vec3 &uv, const vec3 &n, bfpNum etai_over_etat)
{
    bfpNum cos_theta = min(dot(-uv, n), b_1);
    vec3 r_out_perp = etai_over_etat * (uv + cos_theta * n);
    vec3 r_out_parallel = -sqrt(abs(b_1 - r_out_perp.length_squared())) * n;
    return r_out_perp + r_out_parallel;
}

vec3 random_in_unit_disk(bool t)
{
    while (true)
    {
        vec3 p = vec3(random_num(b_1_neg, b_1), random_num(b_1_neg, b_1), b_0);
        if (p.length_squared() >= b_1)
            continue;
        return p;
    }
}

inline void print_vec3(const vec3 &v)
{
    cout << "v[0] : " << bfpNum_to_float(v.e[0]) << endl;
    cout << "v[1] : " << bfpNum_to_float(v.e[1]) << endl;
    cout << "v[2] : " << bfpNum_to_float(v.e[2]) << endl;
}

using point3 = vec3;
using color = vec3;

/* frequently used vec3 */
vec3 vec3_0 = vec3();
point3 point3_0 = vec3();
color color_0 = vec3();

#endif