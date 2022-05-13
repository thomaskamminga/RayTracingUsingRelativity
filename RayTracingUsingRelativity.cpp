#include "rtweekend.h"
#include "ray.h"
#include "vec3.h"
#include "color.h"
#include "hittable_list.h"
#include "sphere.h"
#include "camera.h"
#include "material.h"

#include <iostream>

color ray_color(const ray& r, const hittable_list& world, int depth) {
    hit_record rec;

    //if we've exceeded the ray bounce limit, no more light is gathered.
    if (depth <= 0) {
        return color(0, 0, 0);
    }
    if (world.hit(r, 0.001, infinity, rec)) {
        ray scattered;
        color attenuation;
        if (rec.mat_ptr->scatter(r, rec, attenuation, scattered)) {
            return attenuation * ray_color(scattered, world, depth - 1);
        }
        point3 target = rec.p + rec.normal + random_in_unit_sphere();
        return 0.5 * ray_color(ray(rec.p, target - rec.p), world, depth-1);
    }


    vec3 unit_direction = unit_vector(r.direction());
    auto t = 0.5 * (unit_direction.y() + 1.0);
    return (1.0 - t) * color(1.0, 1.0, 1.0) + t * color(0.5, 0.7, 1.0);

}

hittable_list random_scene() {
    hittable_list world;


    auto checker = make_shared<checker_texture>(color(0.2, 0.3, 0.1), color(0.9, 0.9, 0.9));
    world.add(make_shared<sphere>(point3(0, -1000, 0), 1000, make_shared<lambertian>(checker)));

    for (int a = -11; a < 11; a++) {
        for (int b = -11; b < 11; b++) {
            auto choose_mat = random_double();
            point3 center(a + 0.9 * random_double(), 0.2, b + 0.9 * random_double());

            if ((center - point3(4, 0.2, 0)).length() > 0.9) {
                shared_ptr<material> sphere_material;

                // diffuse
                auto albedo = color::random() * color::random();
                sphere_material = make_shared<lambertian>(albedo);
                world.add(make_shared<sphere>(center, 0.2, sphere_material));
            }
        }
    }

    auto material1 = make_shared<lambertian>(color(0.1, 0.2, 0.5));
    world.add(make_shared<sphere>(point3(0, 1, 0), 1.0, material1));

    auto material2 = make_shared<lambertian>(color(0.4, 0.2, 0.1));
    world.add(make_shared<sphere>(point3(-4, 1, 0), 1.0, material2));

    auto material3 = make_shared<lambertian>(color(0.7, 0.6, 0.5));
    world.add(make_shared<sphere>(point3(4, 1, 0), 1.0, material3));

    return world;
}

hittable_list earth(const point3 lookfrom) {
    hittable_list world;


    auto earth_texture = make_shared<image_texture>("earthmap.jpg");
    auto earth_surface = make_shared<lambertian>(earth_texture);
    auto globe = make_shared<sphere>(point3(0, 0, 0), 2, earth_surface);
    world.add(globe);

    //auto space_texture = make_shared<image_texture>("milkyway_2020_8k.jpg");
    //auto space_surface = make_shared<lambertian>(space_texture);
    //auto space = make_shared<sphere>(lookfrom, 500, space_surface);
    //world.add(space);


    return hittable_list(world);
}

int main()
{
    // Image

    const auto aspect_ratio = 3.0 / 2.0;
    const int image_width = 600;
    const int image_height = static_cast<int>(image_width / aspect_ratio);
    const int samples_per_pixel = 30;
    const int max_depth = 10;

    // Camera
    point3 lookfrom(8, 1, 1);
    point3 lookat(0, 0, 0);
    vec3 vup(0, 1, 0);
    double vof(60);

    camera cam(lookfrom, lookat, vup, vof, aspect_ratio);


    // World

    auto world = earth(lookfrom);

    
    // Render

    std::cout << "P3\n" << image_width << ' ' << image_height << "\n255\n";

    for (int j = image_height-1; j >= 0; --j) {
        std::cerr << "\rScanlines rendered: (" << image_height - j << '/' << image_height << ')' << std::flush;
        for (int i = 0; i < image_width; ++i) {
            color pixel_color(0, 0, 0);
            for (int s = 0; s < samples_per_pixel; s++) {
                auto u = (i + random_double()) / (image_width - 1);
                auto v = (j + random_double()) / (image_height - 1);
                ray r = cam.get_ray(u,v);
                pixel_color += ray_color(r, world, max_depth);
            }
            write_color(std::cout, pixel_color, samples_per_pixel);
        }
    }
    std::cerr << "\nDone.\n";
}


