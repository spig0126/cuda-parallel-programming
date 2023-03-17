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
#include "bvh_bfp.h"
#include "moving_sphere_bfp.h"
#include "sphere_bfp.h"
#include "material_bfp.h"
#include "color_bfp.h"
#include "camera_bfp.h"

#include "mkPpm.h"
#include "mkCuda.h"
#include "mkClockMeasure.h"

#include <ctime>
#include <sys/stat.h>
#include <sys/types.h>
#include <bits/stdc++.h>

#define MAX_SIZE 500

unsigned char *image_array;

// 1. random_scene: Implements the 3D World.
hittable_list random_scene()
{
	hittable_list world;
	int count = 0;

	auto ground_material = make_shared<lambertian>(color(b_0_5, b_0_5, b_0_5));
	world.add(make_shared<sphere>(++count, point3(b_0, b_1000_neg, b_0), b_1000, ground_material));

	auto material1 = make_shared<dielectric>(b_1_5);
	world.add(make_shared<sphere>(++count, point3(b_0, b_1, b_0), b_1, material1));

	auto material2 = make_shared<lambertian>(color(b_0_4, b_0_2, b_0_1));
	world.add(make_shared<sphere>(++count, point3(-b_4, b_1, b_0), b_1, material2));

	auto material3 = make_shared<metal>(color(b_0_7, b_0_6, b_0_5), b_0);
	world.add(make_shared<sphere>(++count, point3(b_4, b_1, b_0), b_1, material3));

	// Constructing BVH
	hittable_list world_bvh;
	world_bvh.add(make_shared<bvh_node>(world, b_0, b_1));

	return world_bvh;
}

// 2. ray_color: calculates color of the current ray intersection point.
color ray_color(const ray &r, const hittable &world, int depth)
{
	hit_record rec;

	if (depth <= 0)
		return color(b_0, b_0, b_0);

	if (world.hit(r, b_0_001, b_infinity, rec)) // hit: 충돌지점을 결정(child ray 의 origin)
	{											// if ray hits object
		ray scattered;
		color attenuation;
		if (rec.mat_ptr->scatter(r, rec, attenuation, scattered)) // scatter: child ray 방향성을 결정
			return attenuation * ray_color(scattered, world, depth - 1);
		return color(b_0, b_0, b_0);
	}

	// if ray doesn;t hit any object: background
	vec3 unit_direction = unit_vector(r.direction());
	bfpNum t = b_0_5 * (unit_direction.y() + b_1);
	return (b_1 - t) * color(b_1, b_1, b_1) + t * color(b_0_5, b_0_7, b_1);
}

color_float ray_color_f(const ray &r, const hittable &world, int depth)
{
	hit_record rec;

	if (depth <= 0)
		return color_float(0, 0, 0);

	if (world.hit(r, b_0_001, b_infinity, rec)) // hit: 충돌지점을 결정(child ray 의 origin)
	{											// if ray hits object
		ray scattered;
		color attenuation;
		if (rec.mat_ptr->scatter(r, rec, attenuation, scattered)) // scatter: child ray 방향성을 결정
			return vec3_bfpNum_to_float(attenuation) * ray_color_f(scattered, world, depth - 1);
		return color_float(0, 0, 0);
	}

	// if ray doesn;t hit any object: background
	vec3_float unit_direction = unit_vector(vec3_bfpNum_to_float(r.direction()));
	float t = 0.5 * (unit_direction.y() + 1.0);
	return (1.0 - t) * color_float(1, 1, 1) + t * color_float(0.5, 0.7, 1.0);
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
	int image_height = static_cast<int>(image_width / aspect_ratio);
	int samples_per_pixel = 10;
	const int max_depth = 50;
	float scale = 1.0 / samples_per_pixel;

	bfpNum _image_height = int_to_bfpNum(image_height);
	bfpNum _image_width = int_to_bfpNum(image_width);
	bfpNum _samples_per_pixel = int_to_bfpNum(samples_per_pixel);

	// World
	ckWolrdBVH->clockResume();
	hittable_list world = random_scene();
	ckWolrdBVH->clockPause();

	// Camera
	point3 lookfrom(b_13, b_2, b_3);
	point3 lookat(b_0, b_0, b_0);
	vec3 vup(b_0, b_1, b_0);
	bfpNum dist_to_focus = b_10;
	bfpNum aperture = b_0_1;
	camera cam(lookfrom, lookat, vup, b_20, float_to_bfpNum(aspect_ratio), aperture, dist_to_focus, b_0, b_1);

	// Rendered Image image_array
	image_array = (unsigned char *)malloc(sizeof(unsigned char) * image_width * image_height * 3);

	// Render
	ckRendering->clockResume();
	float r, g, b;
	for (int j = 0; j < image_height; ++j)
	{
		for (int i = 0; i < image_width; ++i)
		{

			bfpNum _i = int_to_bfpNum(i);
			bfpNum _j = int_to_bfpNum(j);

			int idx = (j * image_width + i) * 3;
			color pixel_color(b_0, b_0, b_0);

			/* float version */
			for (int s = 0; s < samples_per_pixel; ++s)
			{

				bfpNum u = (_i + random_num()) / (_image_width - b_1);
				bfpNum v = ((_image_height - _j - b_1) + random_num()) / (_image_height - b_1);

				// float u = (i + random_float()) / (image_width - 1);
				// float v = ((image_height - j - 1) + random_float()) / (image_height - 1);
				// bfpNum _u = float_to_bfpNum(u);
				// bfpNum _v = float_to_bfpNum(v);

				ray cur_ray = cam.get_ray(u, v);

				pixel_color += ray_color(cur_ray, world, max_depth);
				// pixel_color_f += ray_color_f(cur_ray, world, max_depth);
			}

			r = bfpNum_to_float(pixel_color.x());
			g = bfpNum_to_float(pixel_color.y());
			b = bfpNum_to_float(pixel_color.z());
			// r = pixel_color_f.x();
			// g = pixel_color_f.y();
			// b = pixel_color_f.z();

			// Antialiasing
			// std::cout << "	-> Antialiasing" << endl;

			r = std::sqrt(scale * r);
			g = std::sqrt(scale * g);
			b = std::sqrt(scale * b);

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

	// creating directory
	string directory = "../images/" + to_string(BFP_EXP_BITSIZE) + "_" + to_string(BFP_MANT_BITSIZE);
	if (mkdir(directory.c_str(), 0777) == -1)
		cerr << "Error :  " << strerror(errno) << endl;

	string img_path = directory + "/" + to_string(timeStamp) + "_" + to_string(image_width) + "_" + to_string(samples_per_pixel) + "_" + to_string(max_depth) + "_bfp_img.ppm";
	ppmSave(img_path.c_str(), image_array, image_width, image_height);

	return 0;
}
