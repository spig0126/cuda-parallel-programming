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

unsigned char *array;

// 0. change frequently used int & float values into bfpNum
bfpNum b_0_001 = float_to_bfpNum(0.001);
bfpNum b_0_1 = float_to_bfpNum(0.1);
bfpNum b_0_15 = float_to_bfpNum(0.15);
bfpNum b_0_2 = float_to_bfpNum(0.2);
bfpNum b_0_3 = float_to_bfpNum(0.3);
bfpNum b_0_4 = float_to_bfpNum(0.4);
bfpNum b_0_5 = float_to_bfpNum(0.5);
bfpNum b_0_6 = float_to_bfpNum(0.6);
bfpNum b_0_7 = float_to_bfpNum(0.7);
bfpNum b_0_8 = float_to_bfpNum(0.8);
bfpNum b_0_9 = float_to_bfpNum(0.9);
bfpNum b_0_95 = float_to_bfpNum(0.95);
bfpNum b_1_5 = float_to_bfpNum(1.5);
bfpNum b_1000 = int_to_bfpNum(1000);
bfpNum b_1000_neg = int_to_bfpNum(-1000);
bfpNum b_121 = int_to_bfpNum(121);
bfpNum b_11 = int_to_bfpNum(11);
bfpNum b_3 = int_to_bfpNum(3);
bfpNum b_4 = int_to_bfpNum(4);
bfpNum b_9 = int_to_bfpNum(9);
bfpNum b_10 = int_to_bfpNum(10);
bfpNum b_13 = int_to_bfpNum(13);
bfpNum b_16 = int_to_bfpNum(16);
bfpNum b_20 = int_to_bfpNum(20);

// 1. random_scene: Implements the 3D World.
hittable_list random_scene()
{
	hittable_list world;
	int n = 2;
	int count = 0;

	auto ground_material = make_shared<lambertian>(color(b_0_5, b_0_5, b_0_5));
	world.add(make_shared<sphere>(++count, point3(b_0, b_1000_neg, b_0), b_1000, ground_material));

	for (int a = -n; a < n; a++)
	{
		for (int b = -n; b < n; b++)
		{

			bfpNum a_b = int_to_bfpNum(a);
			bfpNum b_b = int_to_bfpNum(b);

			// Generate constant scene primitives.
			bfpNum choose_mat = (a_b * b_11 + b_b) / b_121;

			// Generate random scene primitives.
			point3 center(a_b, b_0_2, b_b);

			if ((center - point3(b_4, b_0_2, b_0)).length() > b_0_9)
			{

				shared_ptr<material> sphere_material;

				if (choose_mat < b_0_8)
				{

					// diffuse
					color albedo = color::random() * color::random();
					sphere_material = make_shared<lambertian>(albedo);
					point3 center2 = center + vec3(b_0, random_num(b_0, b_0_5), b_0);
					world.add(make_shared<sphere>(++count, center, b_0_2, sphere_material));
				}
				else if (choose_mat < b_0_95)
				{
					// metal
					color albedo = color::random(b_0_5, b_1);
					bfpNum fuzz = random_num(b_0, b_0_5);
					sphere_material = make_shared<metal>(albedo, fuzz);
					world.add(make_shared<sphere>(++count, center, b_0_2, sphere_material));
				}
				else
				{
					// glass
					sphere_material = make_shared<dielectric>(b_1_5);
					world.add(make_shared<sphere>(++count, center, b_0_2, sphere_material));
				}
			}
		}
	}

	auto material1 = make_shared<dielectric>(b_1_5);
	world.add(make_shared<sphere>(++count, point3(b_0, b_1, b_0), b_1, material1));

	auto material2 = make_shared<lambertian>(color(b_0_4, b_0_2, b_0_1));
	world.add(make_shared<sphere>(++count, point3(-b_4, b_1, b_0), b_1, material2));

	auto material3 = make_shared<metal>(color(b_0_7, b_0_6, b_0_5), b_0);
	world.add(make_shared<sphere>(++count, point3(b_4, b_1, b_0), b_1, material3));

	// Constructing BVH
	hittable_list world_bvh;
	world_bvh.add(make_shared<bvh_node>(world, b_0, b_1));
	printf("\n\n================================== BVH CONSTURCTION COMPLETED ==================================\n\n\n");

	return world_bvh;
}

// 2. ray_color: calculates color of the current ray intersection point.
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

	// cout << "Add all sample colors" << endl;
	/* ----------------Add all sample colors------------------------- */
	// for(int i=0; i<ray_colors.size(); i++){
	// 	printf("%f\t%f\t%f\n", ray_colors[i][0], ray_colors[i][1], ray_colors[i][2]);
	// }
	// bfpBlock block = createColorBfpBlock(ray_colors);
	// pixel_color = add_color_bfpBlock(block);
	// r = pixel_color.x();
	// g = pixel_color.y();
	// b = pixel_color.z();
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

	array[idx] = (256 * clamp(r, 0.0, 0.999));
	array[idx + 1] = (256 * clamp(g, 0.0, 0.999));
	array[idx + 2] = (256 * clamp(b, 0.0, 0.999));
}

// 3. main
int main()
{
	// Measure the execution time.
	mkClockMeasure *ckCpu = new mkClockMeasure("CPU CODE");
	ckCpu->clockReset();

	// Image
	auto aspect_ratio = 16.0 / 9.0;
	int image_width = 400; // 400
	int samples_per_pixel = 3;
	const int max_depth = 5;

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
	int image_height = static_cast<int>(image_width / aspect_ratio);
	camera cam(lookfrom, lookat, vup, b_20, float_to_bfpNum(aspect_ratio), aperture, dist_to_focus, b_0, b_1);

	// Rendered Image Array
	unsigned char *bfp_pixel_array = (unsigned char *)malloc(sizeof(unsigned char) * image_width * image_height * 3);

	ckCpu->clockReset();
	ckCpu->clockResume();

	// Render
	for (int j = 0; j < image_height; ++j)
	{
		for (int i = 0; i < image_width; ++i)
		{
			bfp_cal_pixel_color(samples_per_pixel, image_width, image_height, max_depth, world, cam, bfp_pixel_array, i, j);
		}
	}
	cout << "render complete " << endl;

	// RT18 - PRINT PIXEL VALUES OF THE OUTPUT IMAGE:
	// printf("---------------------------------------------\n");
	ckCpu->clockPause();
	ckCpu->clockPrint();

	// for (int i = 0; i < 360000; i += 3)
	// {
	// 	cout << bfp_pixel_array[i] << " " << bfp_pixel_array[i + 1] << " " << bfp_pixel_array[i + 2] << " " << endl;
	// }

	time_t t = time(NULL);
	tm *tPtr = localtime(&t);
	int timeStamp = (((tPtr->tm_year) + 1900) % 100) * 10000 + ((tPtr->tm_mon) + 1) * 100 + (tPtr->tm_mday);
	string bfp_img_path = "../images/" + to_string(BFP_EXP_BITSIZE) + "_" + to_string(BFP_MANT_BITSIZE) + "/" + to_string(timeStamp) + "_" + to_string(image_width) + "_" + to_string(samples_per_pixel) + "_" + to_string(max_depth) + "_bfp_img.ppm";
	ppmSave(bfp_img_path.c_str(), bfp_pixel_array, image_width, image_height);

	return 0;
}
