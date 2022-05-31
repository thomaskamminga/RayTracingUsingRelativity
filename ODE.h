#include <boost/numeric/odeint.hpp>

#ifndef ODE_H
#define ODE_H

using namespace boost::numeric::odeint;

typedef std::vector< double > state_type;


// the system function can be a classical functions
void geodesicODE(state_type& x, state_type& dxdt, double t)
{
	dxdt[0] = x[4];
    dxdt[1] = x[5];
    dxdt[2] = x[6];
    dxdt[3] = x[7];
    dxdt[4] = 0;
    dxdt[5] = 0;
    dxdt[6] = 0;
    dxdt[7] = 0;
}

inline void update_ray(state_type& x, const ray& r) {
    x[0] = r.orig.x();
    x[1] = r.orig.y();
    x[2] = r.orig.z();
    x[3] = 0;
    x[4] = r.dir.x();
    x[5] = r.dir.y();
    x[6] = r.dir.z();
    x[7] = 0;
}

inline void update_ray(ray& r, const state_type& x) {
    point3 origin(x[0], x[1], x[2]);
    point3 direction(x[4], x[5], x[6]);
    r = ray(origin, direction);
}
#endif
