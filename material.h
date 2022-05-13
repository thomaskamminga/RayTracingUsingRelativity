#ifndef MATERIAL_H
#define MATERIAL_H

#include "rtweekend.h"
#include "vec3.h"
#include "ray.h"
#include "hittable.h"
#include "texture.h"

struct hit_record;

class material {
public:
	virtual bool scatter(const ray& r_in, const hit_record& rec, color& attenuation, ray& scattered
	) const = 0;
};


class lambertian : public material {
public:
	lambertian(const color& a) : albedo(make_shared<solid_color>(a)) {}
	lambertian(shared_ptr<texture> a) : albedo(a) {}

	virtual bool scatter(const ray& r_in, const hit_record& rec, color& attenuation, ray& scattered
	) const override {
		auto scatter_direction = rec.normal + random_unit_vector();
		
		// Catch degenerate scatter direction
		if (scatter_direction.near_zero())
			scatter_direction = rec.normal;
		
		scattered = ray(rec.p, scatter_direction);
		attenuation = albedo->value(rec.u, rec.v, rec.p);
		return true;
	}

public:
	shared_ptr<texture> albedo;
};

#endif