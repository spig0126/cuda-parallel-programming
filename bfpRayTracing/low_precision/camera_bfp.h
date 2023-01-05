#ifndef CAMERA_BFP_H
#define CAMERA_BFP_H

#include "ray_bfp.h"

using namespace bfp;

class camera
{
public:
    camera(
        point3 lookfrom,
        point3 lookat,
        vec3 vup,
        bfpNum vfov, // vertical field-of-view in degrees
        bfpNum aspect_ratio,
        bfpNum aperture,
        bfpNum focus_dist,
        bfpNum _time0 = b_0,
        bfpNum _time1 = b_0

    )
    {
        bfpNum theta = degrees_to_radians(vfov);
        bfpNum h = tan(theta / b_2);
        bfpNum viewport_height = b_2 * h;
        bfpNum viewport_width = aspect_ratio * viewport_height;

        vec3 w = unit_vector(lookfrom - lookat);
        vec3 u = unit_vector(cross(vup, w));
        vec3 v = cross(w, u);

        origin = lookfrom;
        horizontal = focus_dist * viewport_width * u;
        vertical = focus_dist * viewport_height * v;
        lower_left_corner = origin - horizontal / b_2 - vertical / b_2 - focus_dist * w;
        lens_radius = aperture / b_2;
        time0 = _time0;
        time1 = _time1;
    }

    ray get_ray(bfpNum s, bfpNum t) const
    {
        vec3 rd = lens_radius * random_in_unit_disk(true);
        vec3 offset = u * rd.x() + v * rd.y();
        return ray(
            origin + offset,
            lower_left_corner + s * horizontal + t * vertical - origin - offset,
            random_num(time0, time1));
    }

public:
    point3 origin;
    point3 lower_left_corner;
    vec3 horizontal;
    vec3 vertical;
    vec3 u, v, w;
    bfpNum lens_radius;
    bfpNum time0, time1; // shutter open/close times
};

#endif