/*
 * ===================================================
 *
 *       Filename:  sphere.h
 *    Description:  Ray Tracing: The Next Week (RTTNW): ~BVH 
 *        Created:  2022/07/13
 * 
 * ===================================================
 */


// Preprocessors

#ifndef SPHERE_H
#define SPHERE_H

#include "hittable.h"
#include "vec3.h"


// Classes

class sphere : public hittable {
public:
    	sphere() {}
        sphere(int idx, point3 cen, float r, shared_ptr<material> m)
            : sphereIdx(idx), center(cen), radius(r), mat_ptr(m) {};
    
    	virtual bool hit(
			const ray& r, float t_min, float t_max, hit_record& rec) const override;

	        virtual bool bounding_box(float time0, float time1, aabb& output_box) const override;

public:
    	point3 center;
    	float radius;
        shared_ptr<material> mat_ptr;
	
	int sphereIdx;  // 구마다 고유번호를 부여!
};


bool sphere::hit(const ray& r, float t_min, float t_max, hit_record& rec) const {
    	vec3 oc = r.origin() - center;
    	auto a = r.direction().length_squared();
    	auto half_b = dot(oc, r.direction());
    	auto c = oc.length_squared() - radius * radius;
    	auto discriminant = half_b * half_b - a * c;

		auto sqrtd = sqrt(discriminant);  // sqrt(b^2-ac)
    	auto root = (-half_b - sqrtd) / a;

		point3 pm = center - vec3(radius, radius, radius);
		point3 pM = center + vec3(radius, radius, radius);


		// If the ray does not hit the sphere, 
    	if (discriminant < 0) {

			// FOR DEBUGGING
			// printf("------------------------------ SPHERE NOT HIT! D/4 = %lf!---------------\n  -%-8s: (%.1lf,%.1lf,%.1lf)\n  -%-8s: (%.1lf, %.1lf, %.1lf) ~ (%lf, %lf, %lf)\n--------------------------------------------------------------------------\n",		     
			// 	discriminant,
			// 	"Ray Position", (r.at(root)).e[0], (r.at(root)).e[1], (r.at(root)).e[2],
	  		// 	"AABB", pm.e[0], pm.e[1], pm.e[2], pM.e[0], pM.e[1], pM.e[2]);
			

			return false;
		}

		// If the ray hits the sphere,
		// FOR DEBUGGING
		// printf("------------------------------ SPHERE HIT!------------------------------\n  -%-8s: (%.1lf,%.1lf,%.1lf)\n  -%-8s: (%.1lf, %.1lf, %.1lf) ~ (%lf, %lf, %lf)\n--------------------------------------------------------------------------\n",
		// "Ray Position", (r.at(root)).e[0], (r.at(root)).e[1], (r.at(root)).e[2],
		// "AABB", pm.e[0], pm.e[1], pm.e[2], pM.e[0], pM.e[1], pM.e[2]);



	//auto sqrtd = sqrt(discriminant);

    	// Find the nearest root that lies in the acceptable range.
    	//auto root = (-half_b - sqrtd) / a;
    	if (root < t_min || t_max < root) {
		root = (-half_b + sqrtd) / a;
		if (root < t_min || t_max < root)
	    		return false;
    	}
    
	rec.t = root;
    	rec.p = r.at(rec.t);
    	vec3 outward_normal = (rec.p - center) / radius;
    	rec.set_face_normal(r, outward_normal);
    	rec.mat_ptr = mat_ptr;


	// printf("------------------------------ HIT RECORD UPDATED ----------------------\n  -%-8s: %lf\n  -%-8s: (%.1lf,%.1lf,%.1lf)\n  -%-8s: (%.1lf, %.1lf, %.1lf) ~ (%lf, %lf, %lf)\n--------------------------------------------------------------------------\n",
	// "Root", rec.t,
	// "Point", (rec.p).e[0], (rec.p).e[1], (rec.p).e[2],
	// "AABB", pm.e[0], pm.e[1], pm.e[2], pM.e[0], pM.e[1], pM.e[2]);
	/*******************************************************/

	/********************** DEBUGGING **********************/
/*
	point3 pm = center - vec3(radius, radius, radius);
	point3 pM = center + vec3(radius, radius, radius);
	printf("***** Sphere %d Hit *****\n", sphereIdx);       
	printf("------------------------------ HIT RECORD ------------------------------\n");
	printf("  -%-8s: %lf\n", "Root", rec.t);
	printf("  -%-8s: (%.1lf,%.1lf,%.1lf)\n", "Point", (rec.p).e[0], (rec.p).e[1], (rec.p).e[2]);
	printf("  -%-8s: (%.1lf, %.1lf, %.1lf) ~ (%lf, %lf, %lf)\n", "AABB", pm.e[0], pm.e[1], pm.e[2], pM.e[0], pM.e[1], pM.e[2]);
	printf("--------------------------------------------------------------------------\n");
*/
	/*******************************************************/




    
	return true;
}


bool sphere::bounding_box(float time0, float time1, aabb& output_box) const {
    	output_box = aabb(
			center - vec3(radius, radius, radius),
			center + vec3(radius, radius, radius));
	
	/********************** DEBUGGING **********************/
	/*
	point3 pm = output_box.minimum, pM = output_box.maximum;
	printf("------------------------------ SPHERE %d INFO ------------------------------\n", sphereIdx);
	printf("  -%-8s: (%.1lf,%.1lf,%.1lf)\n", "Center", center.e[0], center.e[1], center.e[2]);
	printf("  -%-8s: %lf\n", "Radius", radius);
	printf("  -%-8s: (%lf, %lf, %lf) ~ (%lf, %lf, %lf)\n", "AABB", pm.e[0], pm.e[1], pm.e[2], pM.e[0], pM.e[1], pM.e[2]);
	printf("----------------------------------------------------------------------------\n");
	*/
	/*******************************************************/


    	return true;
}


#endif
