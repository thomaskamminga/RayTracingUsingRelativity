#include "mass_object.h"
#include "rtweekend.h"
#include <boost/numeric/odeint.hpp>
#include <cmath>

#ifndef ODE_H
#define ODE_H

using namespace boost::numeric::odeint;

typedef std::vector< double > state_type;

class step {
public:
    double adjusted_step_size(double r) {
        if (use_adjusted_step_size) {
            double temp = clamp(((r) / r_schwarzschild-4),1,infinity) * dt_min;
            return clamp(temp, dt_min, dt_max);
        }
        return dt;
    }
    
    double t;
    double dt;
    double dt_min;
    double dt_max;
    bool use_adjusted_step_size;
    double r_schwarzschild;

};


inline void update_ray(const ray& r_old, ray& r_new, const state_type& x) {

    const point3 new_orig = r_old.orig + r_old.dir;
    const point3 new_dir =  point3(x[0], x[1], x[2]) - new_orig;
    r_new = ray(new_orig, new_dir);
}

inline void update_x(const ray& r, state_type& x) {
    //Update the initial conditions of the ODE to new direction normalize according to the Minkowski metric
    point3 norm_v = (r.dir / r.dir.length())* abs(x[7]);

    x[0] = r.orig[0];
    x[1] = r.orig[1];
    x[2] = r.orig[2];

    x[4] = norm_v[0];
    x[5] = norm_v[1];
    x[6] = norm_v[2];


}

#define Power(x,y) pow(x,y) 
#define Sqrt(x) sqrt(x)


#endif
