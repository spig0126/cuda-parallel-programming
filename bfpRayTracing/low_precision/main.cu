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

#define MAX_SIZE 500

unsigned char *image_array;

// 1. random_scene: Implements the 3D World.
hittable_list random_scene()
{
	// std::cout << "start random_scene()" << endl;

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
	printf("\n================================== BVH CONSTURCTION COMPLETED ==================================\n");

	// std::cout << "end random_scene()" << endl;
	// std::cout << "-----------------------" << endl;

	return world_bvh;
}

// 2. ray_color: calculates color of the current ray intersection point.
color ray_color(const ray &r, const hittable &world, int depth)
{
	// std::cout << "	depth: " << depth << endl;

	hit_record rec;

	if (depth <= 0)
		return color(b_0, b_0, b_0);

	if (world.hit(r, b_0_001, b_infinity, rec))
	{ // if ray hits object
		ray scattered;
		color attenuation;
		if (rec.mat_ptr->scatter(r, rec, attenuation, scattered))
			return attenuation * ray_color(scattered, world, depth - 1);
		return color(b_0, b_0, b_0);
	}

	// if ray doesn;t hit any object: background
	vec3 unit_direction = unit_vector(r.direction());
	bfpNum t = b_0_5 * (unit_direction.y() + b_1);
	return (b_1 - t) * color(b_1, b_1, b_1) + t * color(b_0_5, b_0_7, b_1);
}

color bfp_recur_get_ray_color(const ray &r, const hittable &world, int depth, vector<color> &ray_colors)
{
	// return color(1, 1, 1) when there is color, and color(0,0,0) when max depth has been reached or error
	hit_record rec;

	if (depth <= 0)
		return color(b_0, b_0, b_0);

	if (world.hit(r, b_0_001, b_infinity, rec))
	{ // if ray hits object
		ray scattered;
		color attenuation;
		if (rec.mat_ptr->scatter(r, rec, attenuation, scattered))
		{
			ray_colors.push_back(attenuation);
			bfp_recur_get_ray_color(scattered, world, depth - 1, ray_colors);
			return color(b_1, b_1, b_1);
		}
		else
		{
			return color(b_0, b_0, b_0);
		}
	}
	else
	{ // if ray doesn;t hit any object: background
		vec3 unit_direction = unit_vector(r.direction());
		bfpNum t = b_0_5 * (unit_direction.y() + b_1);
		ray_colors.push_back((b_1 - t) * color(b_1, b_1, b_1) + t * color(b_0_5, b_0_7, b_1));
		return color(b_1, b_1, b_1);
	}
}

color bfp_ray_color(const ray &r, const hittable &world, int depth)
{
	vector<color> ray_colors;

	color res = bfp_recur_get_ray_color(r, world, depth, ray_colors);
	if (isequal(res, color(b_0, b_0, b_0)))
	{
		return res;
	}
	else
	{
		color t = {b_1, b_1, b_1};
		for (int i = 0; i < ray_colors.size(); i++)
		{
			t = t * ray_colors[i];
		}
		return t;
		// return mult_color_bfpBlock(ray_colors);
	}
}

void bfp_cal_pixel_color(int samples_per_pixel, int image_width, int image_height, int max_depth, const hittable &world, camera cam, unsigned char *array, int i, int j)
{
	bfpNum _image_height = int_to_bfpNum(image_height);
	bfpNum _image_width = int_to_bfpNum(image_width);
	bfpNum _i = int_to_bfpNum(i);
	bfpNum _j = int_to_bfpNum(j);

	int idx = (j * image_width + i) * 3;
	color pixel_color(b_0, b_0, b_0);
	float r, g, b;
	vector<color> ray_colors;

	// cout << "Trace sample colors " << endl;
	/* ----------------Trace sample colors ------------------------ */
	for (int s = 0; s < samples_per_pixel; ++s)
	{
		bfpNum u = (_i + random_num()) / (_image_width - b_1);
		bfpNum v = ((_image_height - _j - b_1) + random_num()) / (_image_height - b_1);

		ray cur_ray = cam.get_ray(u, v);

		ray_colors.push_back(bfp_ray_color(cur_ray, world, max_depth)); // add all traced colors to ray_colors
	}

	color c = add_color_bfpBlock(ray_colors);

	r = bfpNum_to_float(c[0]);
	g = bfpNum_to_float(c[1]);
	b = bfpNum_to_float(c[2]);

	// cout << "antialiasing" << endl;
	/* -------------------Antialiasing -----------------------*/
	float scale = 1.0 / samples_per_pixel;
	r = std::sqrt(scale * r);
	g = std::sqrt(scale * g);
	b = std::sqrt(scale * b);

	// cout << r << " " << g << " " << b << endl;

	image_array[idx] = (256 * clamp(r, 0.0, 0.999));
	image_array[idx + 1] = (256 * clamp(g, 0.0, 0.999));
	image_array[idx + 2] = (256 * clamp(b, 0.0, 0.999));
}

// 3. main
int main()
{
	// Measure the execution time.
	mkClockMeasure *ckCpu = new mkClockMeasure("CPU CODE");
	ckCpu->clockReset();

	// Image
	auto aspect_ratio = 16.0 / 9.0;
	int image_width = 1920; // 400
	int image_height = static_cast<int>(image_width / aspect_ratio);
	int samples_per_pixel = 100;
	const int max_depth = 500;
	float scale = 1.0 / samples_per_pixel;

	bfpNum _image_height = int_to_bfpNum(image_height);
	bfpNum _image_width = int_to_bfpNum(image_width);
	bfpNum _samples_per_pixel = int_to_bfpNum(samples_per_pixel);

	ckCpu->clockResume();
	// World
	hittable_list world = random_scene();

	ckCpu->clockPause();
	ckCpu->clockPrint();

	// Camera
	point3 lookfrom(b_13, b_2, b_3);
	point3 lookat(b_0, b_0, b_0);
	vec3 vup(b_0, b_1, b_0);
	bfpNum dist_to_focus = b_10;
	bfpNum aperture = b_0_1;
	camera cam(lookfrom, lookat, vup, b_20, float_to_bfpNum(aspect_ratio), aperture, dist_to_focus, b_0, b_1);

	// Rendered Image image_array
	image_array = (unsigned char *)malloc(sizeof(unsigned char) * image_width * image_height * 3);

	ckCpu->clockReset();
	ckCpu->clockResume();

	// Render
	// std::cout << "start rendering" << endl;
	float r, g, b;
	for (int j = 0; j < image_height; ++j)
	{
		for (int i = 0; i < image_width; ++i)
		{
			// std::cout << j << " " << i << endl;

			bfpNum _i = int_to_bfpNum(i);
			bfpNum _j = int_to_bfpNum(j);

			int idx = (j * image_width + i) * 3;
			color pixel_color(b_0, b_0, b_0);

			/* ----------------Trace sample colors ------------------------ */
			// std::cout << "	-> trace sample colors" << endl;

			for (int s = 0; s < samples_per_pixel; ++s)
			{

				bfpNum u = (_i + random_num()) / (_image_width - b_1);
				bfpNum v = ((_image_height - _j - b_1) + random_num()) / (_image_height - b_1);

				ray cur_ray = cam.get_ray(u, v);

				pixel_color += ray_color(cur_ray, world, max_depth);
			}

			r = bfpNum_to_float(pixel_color.x());
			g = bfpNum_to_float(pixel_color.y());
			b = bfpNum_to_float(pixel_color.z());

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
	// cout << "render complete " << endl;

	// RT18 - PRINT PIXEL VALUES OF THE OUTPUT IMAGE:
	// printf("---------------------------------------------\n");
	ckCpu->clockPause();
	ckCpu->clockPrint();

	// for (int i = 0; i < 360000; i += 3)
	// {
	// 	cout << bfp_pixel_image_array[i] << " " << bfp_pixel_image_array[i + 1] << " " << bfp_pixel_image_array[i + 2] << " " << endl;
	// }

	time_t t = time(NULL);
	tm *tPtr = localtime(&t);
	int timeStamp = (((tPtr->tm_year) + 1900) % 100) * 10000 + ((tPtr->tm_mon) + 1) * 100 + (tPtr->tm_mday);
	string bfp_img_path = "../images/" + to_string(BFP_EXP_BITSIZE) + "_" + to_string(BFP_MANT_BITSIZE) + "/" + to_string(timeStamp) + "_" + to_string(image_width) + "_" + to_string(samples_per_pixel) + "_" + to_string(max_depth) + "_bfp_img.ppm";
	ppmSave(bfp_img_path.c_str(), image_array, image_width, image_height);

	return 0;
}
