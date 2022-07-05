#ifndef HITTABLE_LIST_H
#define HITTABLE_LIST_H

#include "hittable.h"

#include <memory>
#include <vector>

using std::shared_ptr;
using std::make_shared;

class hittable_list : public hittable {
public:
	hittable_list() { background_exits = false; mass_object_exits = false; frame_num = 0;}
	hittable_list(point3 position) { camera_position = position; background_exits = false; mass_object_exits = false; frame_num = 0;}

	void clear() { objects.clear(); }
	void add(shared_ptr<hittable> object) { objects.push_back(object); }
	void add_background(shared_ptr<hittable> object) { background_object = object;  background_exits = true; }
	void add_mass_object(shared_ptr<hittable> object, double mass) { mass_object = object; mass_object_exits = true; }
	void add_frame(point3 position) {
		frame_num++;
		x_cam_vec.push_back(position.x());
		y_cam_vec.push_back(position.y());
		z_cam_vec.push_back(position.z());
		a_vec.push_back(0);
	}
	void add_frame(point3 position, double a) {
		frame_num++;
		x_cam_vec.push_back(position.x());
		y_cam_vec.push_back(position.y());
		z_cam_vec.push_back(position.z());
		a_vec.push_back(a);
	}
	point3 camera_position_frame(int frame) {
		return point3(x_cam_vec[frame], y_cam_vec[frame], z_cam_vec[frame]);
	}
	double rotation_frame(int frame) { return a_vec[frame]; }



	virtual bool hit(
		const ray& r, double t_min, double t_max, hit_record& rec) const override;

public:
	std::vector<shared_ptr<hittable>> objects;
	shared_ptr<hittable> background_object;
	shared_ptr<hittable> mass_object;
	point3 camera_position;
	bool background_exits;
	bool mass_object_exits;
	int frame_num;

private:
	std::vector<double> a_vec;
	std::vector<double> x_cam_vec;
	std::vector<double> y_cam_vec;
	std::vector<double> z_cam_vec;
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
