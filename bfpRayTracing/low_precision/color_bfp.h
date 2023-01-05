#pragma once
#ifndef COLOR_BFP_H
#define COLOR_BFP_H

#include "vec3_bfp.h"

void write_color(std::ostream &out, color pixel_color, int samples_per_pixel)
{
    bfpNum r = pixel_color.x();
    bfpNum g = pixel_color.y();
    bfpNum b = pixel_color.z();

    bfpNum scale = float_to_bfpNum(1.0 / samples_per_pixel);
    r = sqrt(scale * r);
    g = sqrt(scale * g);
    b = sqrt(scale * b);

    out << bfpNum_to_int(b_256 * clamp(r, b_0, b_0_999)) << ' '
        << bfpNum_to_int(b_256 * clamp(g, b_0, b_0_999)) << ' '
        << bfpNum_to_int(b_256 * clamp(b, b_0, b_0_999)) << '\n';
}

void print_color(color c){
    cout << "======== COLOR=======" << endl;
    cout << "r: " << bfpNum_to_float(c[0]) << "\t=>\t" << c[0].sign << " " << c[0].exp << " " << c[0].mant << endl;
    cout << "g: " << bfpNum_to_float(c[1]) << "\t=>\t" << c[1].sign << " " << c[1].exp << " " << c[1].mant << endl;
    cout << "b: " << bfpNum_to_float(c[2]) << "\t=>\t" << c[2].sign << " " << c[2].exp << " " << c[2].mant << endl;
}

color add_color_bfpBlock(vector<color> colors)
{
    bfpNum r, g, b;
    vector<bfpNum> colors_r, colors_g, colors_b;

    // 1. put each r, g, b floats in respective blocks
    for (int i = 0; i < colors.size(); i++)
    {
        colors_r.push_back(colors[i][0]);
        colors_g.push_back(colors[i][1]);
        colors_b.push_back(colors[i][2]);
    }
    bfpBlock r_block = createBfpBlock(colors_r);
    bfpBlock g_block = createBfpBlock(colors_g);
    bfpBlock b_block = createBfpBlock(colors_b);

    // 2. multiply
    r = add_bfpBlock_b(r_block);
    g = add_bfpBlock_b(g_block);
    b = add_bfpBlock_b(b_block);

    return color(r, g, b);
}

color mult_color_bfpBlock(vector<color> colors)
{
    bfpNum r, g, b;
    vector<bfpNum> colors_r, colors_g, colors_b;

    // 1. put each r, g, b floats in respective blocks
    for (int i = 0; i < colors.size(); i++)
    {
        colors_r.push_back(colors[i][0]);
        colors_g.push_back(colors[i][1]);
        colors_b.push_back(colors[i][2]);
    }

    bfpBlock r_block = createBfpBlock(colors_r);
    bfpBlock g_block = createBfpBlock(colors_g);
    bfpBlock b_block = createBfpBlock(colors_b);


    // 2. multiply
    r = mult_bfpBlock_b(r_block);
    g = mult_bfpBlock_b(g_block);
    b = mult_bfpBlock_b(b_block);

    return color(r, g, b);
}



#endif
