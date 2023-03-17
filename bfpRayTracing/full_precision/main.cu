/*
 * ===================================================
 *
 *       Filename:  main.cu
 *    Description:  Ray Tracing In One Weekend (RTIOW): ~BVH
 *        Created:  2022/07/13
 *
 * ===================================================
 */

// Preprocessors

#include "moving_sphere.h"
#include "material.h"
#include "utility.h"
#include "color.h"
#include "hittable_list.h"
#include "sphere.h"
#include "camera.h"
#include "bvh.h"

#include "mkPpm.h"
#include "mkCuda.h"
#include "mkClockMeasure.h"

#include <iostream>
#include <string>

#define MAX_SIZE 500

// using
using namespace std;
unsigned char *image_array;

// Functions

// 1. random_scene: Implements the 3D World.
hittable_list random_scene()
{
	hittable_list world;
	int n = 2;
	int count = 0;

	auto ground_material = make_shared<lambertian>(color(0.5, 0.5, 0.5));
	world.add(make_shared<sphere>(++count, point3(0, -1000, 0), 1000, ground_material));

	auto material1 = make_shared<dielectric>(1.5);
	world.add(make_shared<sphere>(++count, point3(0, 1, 0), 1.0, material1));

	auto material2 = make_shared<lambertian>(color(0.4, 0.2, 0.1));
	world.add(make_shared<sphere>(++count, point3(-4, 1, 0), 1.0, material2));

	auto material3 = make_shared<metal>(color(0.7, 0.6, 0.5), 0.0);
	world.add(make_shared<sphere>(++count, point3(4, 1, 0), 1.0, material3));

	// Constructing BVH
	hittable_list world_bvh;
	world_bvh.add(make_shared<bvh_node>(world, 0, 1));

	return world_bvh;
}

// 2. ray_color: calculates color of the current ray intersection point.
color ray_color(const ray &r, const hittable &world, int depth)
{

	hit_record rec;

	// Limit the number of child ray.
	if (depth <= 0)
		return color(0, 0, 0); // If the ray hits objects more than 'depth' times, consider that no light approaches the current point.

	// If the ray hits an object: Hittable Object
	if (world.hit(r, 0.001, infinity, rec))
	{

		ray scattered;
		color attenuation;
		if (rec.mat_ptr->scatter(r, rec, attenuation, scattered))
			return attenuation * ray_color(scattered, world, depth - 1);
		return color(0, 0, 0);
	}

	// If the ray hits no object: Background
	vec3 unit_direction = unit_vector(r.direction());
	auto t = 0.5 * (unit_direction.y() + 1.0);
	return (1.0 - t) * color(1.0, 1.0, 1.0) + t * color(0.5, 0.7, 1.0);
}

// 3. main
int main()
{
	// Measure the execution time.
	mkClockMeasure *ckCpu = new mkClockMeasure("TOTAL TIME");
	mkClockMeasure *ckWolrdBVH = new mkClockMeasure("CONSTRUCT WORLD & BVH");
	mkClockMeasure *ckRendering = new mkClockMeasure("RENDERING");
	ckCpu->clockReset();
	ckWolrdBVH->clockReset();
	ckRendering->clockReset();

	ckCpu->clockResume();
	// Image
	auto aspect_ratio = 16.0 / 9.0;
	int image_width = 400; // 400
	int samples_per_pixel = 10;
	const int max_depth = 50;

	// World
	ckWolrdBVH->clockResume();
	hittable_list world = random_scene();
	ckWolrdBVH->clockPause();

	// Camera
	point3 lookfrom(13, 2, 3);
	point3 lookat(0, 0, 0);
	vec3 vup(0, 1, 0);
	auto dist_to_focus = 10.0;
	auto aperture = 0.1;
	int image_height = static_cast<int>(image_width / aspect_ratio);
	camera cam(lookfrom, lookat, vup, 20, aspect_ratio, aperture, dist_to_focus, 0.0, 1.0);

	// Rendered Image Array
	image_array = (unsigned char *)malloc(sizeof(unsigned char) * image_width * image_height * 3);

	// Render
	ckRendering->clockResume();
	float r, g, b;

	for (int j = 0; j < image_height; ++j)
	{
		for (int i = 0; i < image_width; ++i)
		{
			int idx = (j * image_width + i) * 3;
			color pixel_color(0, 0, 0);

			for (int s = 0; s < samples_per_pixel; ++s)
			{
				auto u = (i + random_float()) / (image_width - 1);
				auto v = ((image_height - j - 1) + random_float()) / (image_height - 1);

				ray cur_ray = cam.get_ray(u, v);
				
				pixel_color += ray_color(cur_ray, world, max_depth);

				r = pixel_color.x();
				g = pixel_color.y();
				b = pixel_color.z();

				// Antialiasing
				float scale = 1.0 / samples_per_pixel;
				r = sqrt(scale * r);
				g = sqrt(scale * g);
				b = sqrt(scale * b);
			}

			image_array[idx] = (256 * clamp(r, 0.0, 0.999));
			image_array[idx + 1] = (256 * clamp(g, 0.0, 0.999));
			image_array[idx + 2] = (256 * clamp(b, 0.0, 0.999));
		}
	}
	ckRendering->clockPause();
	ckCpu->clockPause();

	/* Print clocks */
	ckCpu->clockPrint();
	ckWolrdBVH->clockPrint();
	ckRendering->clockPrint();

	/* Save image */
	time_t t = time(NULL);
	tm *tPtr = localtime(&t);
	int timeStamp = (((tPtr->tm_year) + 1900) % 100) * 10000 + ((tPtr->tm_mon) + 1) * 100 + (tPtr->tm_mday);
	string img_path = "../images/FP32/" + to_string(timeStamp) + "_" + to_string(image_width) + "_" + to_string(samples_per_pixel) + "_" + to_string(max_depth) + "_img.ppm";

	ppmSave(img_path.c_str(), image_array, image_width, image_height);
	return 0;
}