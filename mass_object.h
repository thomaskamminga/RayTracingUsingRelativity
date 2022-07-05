#ifndef MASS_OBJECT_H
#define MASS_OBJECT_H

#include "vec3.h"

class mass_object {
public:
	mass_object(point3 location, double mass) {
		loc = location;
		m = mass;
	}
	
	point3 location() { return loc; }
	double mass() { return m;  }

private:
	point3 loc;
	double m;
};


#endif