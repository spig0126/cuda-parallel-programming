/*
 * ===================================================
 *      
 *       Filename:  camera.h
 *    Description:  Ray Tracing: The Next Week (RTTNW): ~BVH 
 *        Created:  2022/07/13
 * 
 * ===================================================
 */


// Preprocessors

#ifndef CAMERA_H
#define CAMERA_H
#include "utility.h"


// Classes

class camera {
public:
	camera(
    		point3 lookfrom,
    		point3 lookat,
    		vec3   vup,	
		float vfov, // vertical field-of-view in degrees
		float aspect_ratio,
    		float aperture,
		float focus_dist,
		float _time0 = 0,
    		float _time1 = 0

	      ) {
		auto theta = degrees_to_radians(vfov);
		auto h = tan(theta/2);
    		auto viewport_height = 2.0 * h;
    		auto viewport_width = aspect_ratio * viewport_height;		

    		auto w = unit_vector(lookfrom - lookat);
    		auto u = unit_vector(cross(vup, w));
    		auto v = cross(w, u);

    		origin = lookfrom;
		horizontal = focus_dist * viewport_width * u;
    		vertical = focus_dist * viewport_height * v;
		lower_left_corner = origin - horizontal/2 - vertical/2 - focus_dist*w;
	    	lens_radius = aperture / 2;
	       	time0 = _time0;
    		time1 = _time1;

	}

        ray get_ray(float s, float t) const {
	  	vec3 rd = lens_radius * random_in_unit_disk();
	      	vec3 offset = u * rd.x() + v * rd.y();
    		return ray(
		    	origin + offset,
			lower_left_corner + s*horizontal + t*vertical - origin - offset,
			random_float(time0, time1)
		);

	}

private:
	point3 origin;
    	point3 lower_left_corner;
    	vec3 horizontal;
    	vec3 vertical;
	vec3 u, v, w;
        float lens_radius;
        float time0, time1;  // shutter open/close times

};


#endif

