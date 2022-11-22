/*
 * ===================================================
 *
 *       Filename:  hittable.h
 *    Description:  Ray Tracing: The Next Week (RTTNW): ~BVH 
 *        Created:  2022/07/13
 * 
 * ===================================================
 */


// Preprocessors

#ifndef HITTABLE_H
#define HITTABLE_H

#include "aabb.h"
#include "ray.h"
#include "utility.h"


class material;

// Classes
struct hit_record {
    	point3 p;
    	vec3 normal;
	shared_ptr<material> mat_ptr;
    	float t;
       	bool front_face;

    	inline void set_face_normal(const ray& r, const vec3& outward_normal) {
		front_face = dot(r.direction(), outward_normal) < 0;
		normal = front_face ? outward_normal : -outward_normal;
    	}
};


class hittable {
	public:
		virtual bool hit(const ray& r, float t_min, float t_max, hit_record& rec) const = 0;
		virtual bool bounding_box(float time0, float time1, aabb& output_box) const = 0;

};


#endif
