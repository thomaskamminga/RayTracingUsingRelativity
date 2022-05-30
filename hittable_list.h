#ifndef HITTABLE_LIST_H
#define HITTABLE_LIST_H

#include "hittable.h"

#include <memory>
#include <vector>

using std::shared_ptr;
using std::make_shared;

class hittable_list : public hittable {
	public:
		hittable_list() { background_exits = false; }
		hittable_list(point3 position) { camera_position = position; background_exits = false;}

		void clear() { objects.clear(); }
		void add(shared_ptr<hittable> object) { objects.push_back(object); }
		void add_background(shared_ptr<hittable> object) { background_object = object;  background_exits = true;}


		virtual bool hit(
			const ray& r, double t_min, double t_max, hit_record& rec) const override;

	public:
		std::vector<shared_ptr<hittable>> objects;
		shared_ptr<hittable> background_object;
		point3 camera_position;
		bool background_exits;
};



bool hittable_list::hit(const ray& r, double t_min, double t_max, hit_record& rec) const {
	hit_record temp_rec;
	bool hit_anything = false;
	auto closest_so_far = t_max;

	for (const auto& object : objects) {
		if (object->hit(r, t_min, closest_so_far, temp_rec)) {
			hit_anything = true;
			closest_so_far = temp_rec.t;
			rec = temp_rec;
		}
	}
	return hit_anything;
}




#endif // !HITTABLE_LIST_H
