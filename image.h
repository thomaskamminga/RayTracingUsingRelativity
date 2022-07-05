#ifndef COLOR_H
#define COLOR_H

#include "vec3.h"

#include <iostream>

void write_color(unsigned char image_data[], color pixel_color, int samples_per_pixel, int i, int j, int image_width, int image_height, int channel_num) {
    auto r = pixel_color.x();
    auto g = pixel_color.y();
    auto b = pixel_color.z();
    
    // Divide the color by the number of samples and gamma-correct for gamma=2.0.
    auto scale = 1.0 / samples_per_pixel;
    r = sqrt(scale * r);
    g = sqrt(scale * g);
    b = sqrt(scale * b);

    // Write the translated [0,255] value of each color component.
    const int position = channel_num * ((image_height - j - 1) * image_width + i );
    image_data[position + 0] = static_cast<unsigned char>(256 * clamp(r, 0.0, 0.999));
    image_data[position + 1] = static_cast<unsigned char>(256 * clamp(g, 0.0, 0.999));
    image_data[position + 2] = static_cast<unsigned char>(256 * clamp(b, 0.0, 0.999));
}

#endif
