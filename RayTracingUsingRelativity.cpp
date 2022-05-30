#include "rtweekend.h"
#include "ray.h"
#include "vec3.h"
#include "color.h"
#include "hittable_list.h"
#include "sphere.h"
#include "camera.h"
#include "material.h"
#include "ODE.h"

#include <boost/numeric/odeint.hpp>
#include <iostream>





ray next_ray_on_geodesic(const ray& r, const double step_size) {
    ray next_ray;

    next_ray.orig = r.orig + step_size * r.dir;
    next_ray.dir = r.dir;
    return next_ray;
}

color ray_color(const ray& r, const hittable_list& world, const double step_size, const double background_radius, int bounce_depth, int step_depth) {
    hit_record rec;

    //if we've exceeded the ray bounce limit, no more light is gathered.
    if (bounce_depth <= 0) {
        return color(0, 0, 0);
    }

    //if we exceed search length return pink for debugging (avoids infinte loop)
    if (step_depth <= 0) {
        return color(1, 0, 1);
    }

    //ray hits object
    if (world.hit(r, 0.001, step_size, rec)) {
        ray scattered;
        color attenuation;

        if (rec.mat_ptr->scatter(r, rec, attenuation, scattered)) {
            return 5 * attenuation * ray_color(scattered, world, step_size, background_radius, bounce_depth - 1, step_depth);
        }
        return color(0, 0, 0);
    }

    
    //ray does not hit object and is out of bound
    if ((r.at(step_size) - world.camera_position).length() > background_radius) {
        if (world.background_exits) {
            if (world.background_object->hit(ray(point3(0, 0, 0), r.dir), 0.001, infinity, rec)) {
                ray scattered;
                color attenuation;

                if (rec.mat_ptr->scatter(r, rec, attenuation, scattered)) {
                    return attenuation;
                }
            }
        }
        // no background object 
        vec3 unit_direction = unit_vector(r.direction());
        auto t = 0.5 * (unit_direction.y() + 1.0);
        return (1.0 - t) * color(1.0, 1.0, 1.0) + t * color(0.5, 0.7, 1.0);
    }

    //ray does not hit object and is in bounds
    ray next_ray = next_ray_on_geodesic(r, step_size);
    return ray_color(next_ray, world, step_size, background_radius, bounce_depth, step_depth-1);


}

hittable_list random_scene(const point3 lookfrom) {
    hittable_list world(lookfrom);


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
    hittable_list world(lookfrom);


    auto earth_texture = make_shared<image_texture>("earthmap.jpg");
    auto earth_surface = make_shared<lambertian>(earth_texture);
    auto globe = make_shared<sphere>(point3(0, 0, 0), 2, earth_surface);
    world.add(globe);

    auto space_texture = make_shared<image_texture>("milkyway_2020_8k.jpg");
    auto space_surface = make_shared<lambertian>(space_texture);
    auto space = make_shared<sphere>(lookfrom, 500, space_surface);
    world.add_background(space);


    return hittable_list(world);
}


int main()
{
    // Image

    const auto aspect_ratio = 16.0 / 9.0;
    const int image_width = 1920;
    const int image_height = static_cast<int>(image_width / aspect_ratio);
    const int samples_per_pixel = 200;
    const int max_bounce_depth = 10;

    // Geodesic

    const auto background_radius = 10.0;
    const auto max_geodesic_length = 70.0;
    const int geodesic_depth = 10;
    const auto geodesic_step_size = static_cast<double>(max_geodesic_length / geodesic_depth);

    // Camera

    point3 lookfrom(0, 0, -8);
    point3 lookat(0, 0, 0);
    vec3 vup(0, 1, 0);
    double vof(100);

    camera cam(lookfrom, lookat, vup, vof, aspect_ratio);


    // World

    auto world = earth(lookfrom);


    // OdeInt
    
    //state_type x(3);
    //x[0] = x[1] = x[2] = 10.0;
    //const double dt = 0.01;
    //runge_kutta4< state_type > rk4;
    //double t = 0.0;
    //for (size_t i = 0; i < 1000; ++i, t += dt)
    //{
    //    rk4.do_step(lorenz, x, t, dt);
    //    std::cerr << x[0] << ' ' << x[1] << ' ' << x[3] << ' ' << '\n';
    //    std::cin.ignore();
    //}
    //


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
                pixel_color += ray_color(r, world, geodesic_step_size, background_radius, max_bounce_depth, geodesic_depth);
            }
            write_color(std::cout, pixel_color, samples_per_pixel);
        }
    }
    std::cerr << "\nDone.\n";
}


