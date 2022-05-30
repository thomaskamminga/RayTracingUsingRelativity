#include <boost/numeric/odeint.hpp>

#ifndef ODE_H
#define ODE_H

using namespace boost::numeric::odeint;

typedef std::vector< double > state_type;

const double sigma = 10.0;
const double R = 28.0;
const double b = 8.0 / 3.0;

// the system function can be a classical functions
void lorenz(state_type& x, state_type& dxdt, double t)
{
    dxdt[0] = sigma * (x[1] - x[0]);
    dxdt[1] = R * x[0] - x[1] - x[0] * x[2];
    dxdt[2] = x[0] * x[1] - b * x[2];
}

// the system function can also be a functor
class lorenz_class
{
public:
    void operator()(state_type& x, state_type& dxdt, double t)
    {
        dxdt[0] = sigma * (x[1] - x[0]);
        dxdt[1] = R * x[0] - x[1] - x[0] * x[2];
        dxdt[2] = x[0] * x[1] - b * x[2];
    }
};


#endif
