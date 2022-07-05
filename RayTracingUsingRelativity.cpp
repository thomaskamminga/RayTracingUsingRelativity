#include "rtweekend.h"
#include "ray.h"
#include "vec3.h"
#include "image.h"
#include "hittable_list.h"
#include "sphere.h"
#include "camera.h"
#include "material.h"
#include "ODE.h"
#include "mass_object.h"
#include "aarect.h"

#include <boost/numeric/odeint.hpp>
#include <iostream>
#include <string>




// ODE variables

const int metric_case = 1;
const auto M = 1.0;
const auto G = 1.0;
double a;
const auto r_schwarzschild = 2 * G * M;
auto star = mass_object(point3(0, 0, 0), M);

void geodesicODE(const state_type& x, state_type& dxdt, const double t)
{
    double r;
    double r_X, r_Y, r_Z;
    double f_xyz;

    switch (metric_case) 
    {
    case 0: // No curvature
        dxdt[0] = x[4];
        dxdt[1] = x[5];
        dxdt[2] = x[6];
        dxdt[3] = x[7];
        dxdt[4] = 0;
        dxdt[5] = 0;
        dxdt[6] = 0;
        dxdt[7] = 0;

        break;
                 
    case 1: // Kerr-Schild no rotation
        r = (point3(x[0], x[1], x[2]) - star.location()).length();

        f_xyz = (G * M * ((2 * G * M - r) * Power(r, 2) * Power(x[7], 2) + (-2 * Power(r, 3) + (2 * G * M + 3 * r) * Power(x[0], 2)) * Power(x[4], 2) + (2 * G * M * Power(x[1], 2) + r * (-2 * Power(x[0], 2) + Power(x[1], 2) - 2 * Power(x[2], 2))) * Power(x[5], 2) + 2 * (2 * G * M + 3 * r) * x[1] * x[2] * x[5] * x[6] + ((2 * G * M + r) * Power(x[2], 2) - 2 * r * (Power(r, 2) - Power(x[2], 2))) * Power(x[6], 2) + 2 * (2 * G * M + 3 * r) * x[0] * x[4] * (x[1] * x[5] + x[2] * x[6]) + 4 * G * M * r * x[7] * (x[0] * x[4] + x[1] * x[5] + x[2] * x[6]))) / Power(r, 6);
        dxdt[0] = x[4];
        dxdt[1] = x[5];
        dxdt[2] = x[6];
        dxdt[3] = x[7];
        dxdt[4] = x[0] * f_xyz;
        dxdt[5] = x[1] * f_xyz;
        dxdt[6] = x[2] * f_xyz;
        dxdt[7] = (2 * G * M * (-(G * M * Power(r, 2) * Power(x[7], 2)) + (-(G * M * Power(x[0], 2)) + r * (Power(r, 2) - 2 * Power(x[0], 2))) * Power(x[4], 2) + (-(G * M * Power(x[1], 2)) + r * (Power(r, 2) - 2 * Power(x[1], 2))) * Power(x[5], 2) - 2 * (G * M + 2 * r) * x[1] * x[2] * x[5] * x[6] + (-((G * M + r) * Power(x[2], 2)) + r * (Power(r, 2) - Power(x[2], 2))) * Power(x[6], 2) - 2 * (G * M + 2 * r) * x[0] * x[4] * (x[1] * x[5] + x[2] * x[6]) - (2 * G * M * r + Power(r, 2)) * x[7] * (x[0] * x[4] + x[1] * x[5] + x[2] * x[6]))) / Power(r, 5);

        break;
    case 2:  // Kerr-Schild with rotation
        r =   Sqrt(-Power(a, 2) + Power(x[0], 2) + Power(x[1], 2) + Power(x[2], 2) + \
                Sqrt(4 * Power(a, 2) * Power(x[2], 2) + Power(-Power(a, 2) + Power(x[0], 2) + \
                Power(x[1], 2) + Power(x[2], 2), 2))) / Sqrt(2);
        r_X = (x[0]*Sqrt(-Power(a,2) + Power(x[0],2) + Power(x[1],2) + \
                Power(x[2],2) + Sqrt(4*Power(a,2)*Power(x[2],2) + Power(-Power(a,2) + \
                Power(x[0],2) + Power(x[1],2) + \
                Power(x[2],2),2))))/(Sqrt(2)*Sqrt(4*Power(a,2)*Power(x[2],2) + \
                Power(-Power(a,2) + Power(x[0],2) + Power(x[1],2) + \
                Power(x[2],2),2)));
        r_Y = (x[1]*Sqrt(-Power(a,2) + Power(x[0],2) + Power(x[1],2) + \
                Power(x[2],2) + Sqrt(4*Power(a,2)*Power(x[2],2) + Power(-Power(a,2) + \
                Power(x[0],2) + Power(x[1],2) + \
                Power(x[2],2),2))))/(Sqrt(2)*Sqrt(4*Power(a,2)*Power(x[2],2) + \
                Power(-Power(a,2) + Power(x[0],2) + Power(x[1],2) + \
                Power(x[2],2),2)));
        r_Z = (x[2]*(1 + (Power(a,2) + Power(x[0],2) + Power(x[1],2) + \
                Power(x[2],2))/Sqrt(4*Power(a,2)*Power(x[2],2) + Power(-Power(a,2) + \
                Power(x[0],2) + Power(x[1],2) + \
                Power(x[2],2),2))))/(Sqrt(2)*Sqrt(-Power(a,2) + Power(x[0],2) + \
                Power(x[1],2) + Power(x[2],2) + Sqrt(4*Power(a,2)*Power(x[2],2) + \
                Power(-Power(a,2) + Power(x[0],2) + Power(x[1],2) + \
                Power(x[2],2),2))));

        dxdt[0] = x[4];
        dxdt[1] = x[5];
        dxdt[2] = x[6];
        dxdt[3] = x[7];
        dxdt[4] = (G*M*(-(Power(x[7],2)*Power(Power(a,2) + \
Power(r,2),3)*(-2*G*M*x[2]*Power(r,4)*(a*x[1] + x[0]*r)*(Power(a,2) + \
Power(r,2))*(2*Power(a,2)*x[2]*r - 3*Power(a,2)*Power(x[2],2)*r_Z + \
Power(r,4)*r_Z) - 2*G*M*Power(r,3)*(a*x[1] + x[0]*r)*(-(a*x[0]) + \
x[1]*r)*(-3*Power(a,2)*Power(x[2],2)*Power(r,2) + Power(r,6))*r_Y - \
(-3*Power(a,2)*Power(x[2],2)*Power(r,2) + \
Power(r,6))*(2*G*M*Power(r,3)*Power(a*x[1] + x[0]*r,2) - \
Power(Power(a,2) + Power(r,2),2)*(Power(a,2)*Power(x[2],2) + \
Power(r,4)))*r_X)) - 2*x[7]*x[6]*r*(Power(a,2) + \
Power(r,2))*(-2*G*M*Power(r,4)*(a*x[1] + x[0]*r)*Power(Power(a,2) + \
Power(r,2),3)*(2*Power(a,2)*x[2]*r - 3*Power(a,2)*Power(x[2],2)*r_Z + \
Power(r,4)*r_Z) + 2*G*M*Power(r,3)*(a*x[1] + x[0]*r)*(a*x[0] - \
x[1]*r)*(2*Power(r,3)*(a*x[0] - x[1]*r)*(Power(a,2)*Power(x[2],2) + \
Power(r,4))*r_Z + x[1]*Power(r,2)*(Power(a,2) + \
Power(r,2))*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_Z + \
3*r*(-(a*x[0]) + x[1]*r)*(Power(a,2) + \
Power(r,2))*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_Z + \
Power(r,2)*(a*x[0] - x[1]*r)*(Power(a,2) + \
Power(r,2))*(2*Power(a,2)*x[2] + 4*Power(r,3)*r_Z) + \
4*x[2]*Power(r,4)*Power(Power(a,2) + Power(r,2),2)*r_Y - \
2*x[2]*Power(Power(a,2) + Power(r,2),2)*(Power(a,2)*Power(x[2],2) + \
Power(r,4))*r_Y) + (2*G*M*Power(r,3)*Power(a*x[1] + x[0]*r,2) - \
Power(Power(a,2) + Power(r,2),2)*(Power(a,2)*Power(x[2],2) + \
Power(r,4)))*(2*Power(r,3)*(a*x[1] + \
x[0]*r)*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_Z - \
x[0]*Power(r,2)*(Power(a,2) + Power(r,2))*(Power(a,2)*Power(x[2],2) + \
Power(r,4))*r_Z - 3*r*(a*x[1] + x[0]*r)*(Power(a,2) + \
Power(r,2))*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_Z + \
2*Power(r,2)*(a*x[1] + x[0]*r)*(Power(a,2) + \
Power(r,2))*(Power(a,2)*x[2] + 2*Power(r,3)*r_Z) - \
4*x[2]*Power(r,4)*Power(Power(a,2) + Power(r,2),2)*r_X + \
2*x[2]*Power(Power(a,2) + Power(r,2),2)*(Power(a,2)*Power(x[2],2) + \
Power(r,4))*r_X)) + Power(x[6],2)*(Power(a,2) + \
Power(r,2))*(-4*G*M*Power(r,4)*(a*x[1] + x[0]*r)*Power(Power(a,2) + \
Power(r,2),3)*(-(Power(a,2)*Power(x[2],2)) + Power(r,4))*(r - \
2*x[2]*r_Z) + 2*G*M*Power(x[2],2)*Power(r,2)*(a*x[1] + \
x[0]*r)*Power(Power(a,2) + Power(r,2),3)*(2*Power(r,5) + \
Power(a,2)*Power(x[2],3)*r_Z - 3*x[2]*Power(r,4)*r_Z) + \
2*G*M*Power(r,3)*(a*x[1] + x[0]*r)*(a*x[0] - \
x[1]*r)*(2*Power(a,4)*x[1]*Power(x[2],2)*Power(r,3) - \
2*x[1]*Power(r,9) + 4*Power(a,5)*x[0]*Power(x[2],3)*r*r_Z + \
2*Power(a,2)*x[2]*Power(r,5)*(x[1]*x[2] - 2*a*x[0]*r_Z) - \
2*a*Power(r,7)*(a*x[1] + 4*x[0]*x[2]*r_Z) + \
Power(a,6)*Power(x[2],4)*r_Y + 2*Power(a,2)*Power(r,6)*(a*x[0] + \
x[1]*x[2]*r_Z - 3*Power(x[2],2)*r_Y) + Power(r,8)*(2*a*x[0] + \
6*x[1]*x[2]*r_Z - 3*Power(x[2],2)*r_Y) - \
2*Power(a,4)*Power(x[2],2)*Power(r,2)*(a*x[0] + 3*x[1]*x[2]*r_Z - \
Power(x[2],2)*r_Y) + Power(a,2)*Power(x[2],2)*Power(r,4)*(-2*a*x[0] - \
2*x[1]*x[2]*r_Z + (-3*Power(a,2) + Power(x[2],2))*r_Y)) + \
(2*G*M*Power(r,3)*Power(a*x[1] + x[0]*r,2) - Power(Power(a,2) + \
Power(r,2),2)*(Power(a,2)*Power(x[2],2) + \
Power(r,4)))*(2*Power(r,2)*(a*x[1] + x[0]*r)*(Power(a,2) + \
Power(r,2))*(Power(a,2)*Power(x[2],2) + Power(r,4)) - \
4*x[2]*Power(r,3)*(a*x[1] + x[0]*r)*(Power(a,2)*Power(x[2],2) + \
Power(r,4))*r_Z + 2*x[0]*x[2]*Power(r,2)*(Power(a,2) + \
Power(r,2))*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_Z + \
4*x[2]*r*(a*x[1] + x[0]*r)*(Power(a,2) + \
Power(r,2))*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_Z - \
4*x[2]*Power(r,2)*(a*x[1] + x[0]*r)*(Power(a,2) + \
Power(r,2))*(Power(a,2)*x[2] + 2*Power(r,3)*r_Z) + \
4*Power(x[2],2)*Power(r,4)*Power(Power(a,2) + Power(r,2),2)*r_X - \
Power(x[2],2)*Power(Power(a,2) + \
Power(r,2),2)*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_X)) + \
4*G*M*x[7]*x[4]*Power(r,3)*(a*x[1] + x[0]*r)*(Power(a,2) + \
Power(r,2))*(Power(Power(a,2) + \
Power(r,2),3)*(-3*Power(a,2)*Power(x[2],2)*Power(r,2) + \
Power(r,6))*r_X + x[2]*(Power(a,2) + \
Power(r,2))*(2*Power(r,3)*(a*x[1] + x[0]*r)*(Power(a,2)*Power(x[2],2) \
+ Power(r,4))*r_Z - x[0]*Power(r,2)*(Power(a,2) + \
Power(r,2))*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_Z - 3*r*(a*x[1] \
+ x[0]*r)*(Power(a,2) + Power(r,2))*(Power(a,2)*Power(x[2],2) + \
Power(r,4))*r_Z + 2*Power(r,2)*(a*x[1] + x[0]*r)*(Power(a,2) + \
Power(r,2))*(Power(a,2)*x[2] + 2*Power(r,3)*r_Z) - \
4*x[2]*Power(r,4)*Power(Power(a,2) + Power(r,2),2)*r_X + \
2*x[2]*Power(Power(a,2) + Power(r,2),2)*(Power(a,2)*Power(x[2],2) + \
Power(r,4))*r_X) + Power(r,2)*(a*x[0] - \
x[1]*r)*(2*Power(a,3)*Power(r,5) + \
3*Power(a,5)*Power(x[2],2)*(x[1]*r_Y + x[0]*r_X) + \
Power(a,3)*Power(x[2],2)*Power(r,2)*(x[1]*r_Y + x[0]*r_X) - \
Power(a,3)*Power(r,4)*(x[1]*r_Y + x[0]*r_X) - \
3*a*Power(r,6)*(x[1]*r_Y + x[0]*r_X) + \
2*Power(a,4)*Power(x[2],2)*r*(a + 2*x[0]*r_Y - 2*x[1]*r_X) + \
2*Power(a,2)*Power(x[2],2)*Power(r,3)*(a + x[0]*r_Y - x[1]*r_X) + \
2*Power(r,7)*(a - x[0]*r_Y + x[1]*r_X))) + \
2*x[7]*x[5]*Power(r,2)*(Power(a,2) + Power(r,2))*(2*G*M*r*(a*x[1] + \
x[0]*r)*Power(Power(a,2) + \
Power(r,2),3)*(-3*Power(a,2)*Power(x[2],2)*Power(r,2) + \
Power(r,6))*r_Y - 2*G*M*x[2]*r*(a*x[1] + x[0]*r)*(Power(a,2) + \
Power(r,2))*(2*Power(r,3)*(a*x[0] - x[1]*r)*(Power(a,2)*Power(x[2],2) \
+ Power(r,4))*r_Z + x[1]*Power(r,2)*(Power(a,2) + \
Power(r,2))*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_Z + \
3*r*(-(a*x[0]) + x[1]*r)*(Power(a,2) + \
Power(r,2))*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_Z + \
Power(r,2)*(a*x[0] - x[1]*r)*(Power(a,2) + \
Power(r,2))*(2*Power(a,2)*x[2] + 4*Power(r,3)*r_Z) + \
4*x[2]*Power(r,4)*Power(Power(a,2) + Power(r,2),2)*r_Y - \
2*x[2]*Power(Power(a,2) + Power(r,2),2)*(Power(a,2)*Power(x[2],2) + \
Power(r,4))*r_Y) + (2*G*M*Power(r,3)*Power(a*x[1] + x[0]*r,2) - \
Power(Power(a,2) + Power(r,2),2)*(Power(a,2)*Power(x[2],2) + \
Power(r,4)))*(2*Power(a,3)*Power(r,5) + \
3*Power(a,5)*Power(x[2],2)*(x[1]*r_Y + x[0]*r_X) + \
Power(a,3)*Power(x[2],2)*Power(r,2)*(x[1]*r_Y + x[0]*r_X) - \
Power(a,3)*Power(r,4)*(x[1]*r_Y + x[0]*r_X) - \
3*a*Power(r,6)*(x[1]*r_Y + x[0]*r_X) + \
2*Power(a,4)*Power(x[2],2)*r*(a + 2*x[0]*r_Y - 2*x[1]*r_X) + \
2*Power(a,2)*Power(x[2],2)*Power(r,3)*(a + x[0]*r_Y - x[1]*r_X) + \
2*Power(r,7)*(a - x[0]*r_Y + x[1]*r_X))) + \
Power(x[4],2)*Power(r,2)*(a*x[1] + \
x[0]*r)*(-((2*G*M*Power(r,3)*Power(a*x[1] + x[0]*r,2) - \
Power(Power(a,2) + Power(r,2),2)*(Power(a,2)*Power(x[2],2) + \
Power(r,4)))*(4*Power(r,4)*(a*x[1] + x[0]*r)*(Power(a,2) + \
Power(r,2))*r_X + 4*Power(r,2)*(a*x[1] + \
x[0]*r)*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_X - 3*(a*x[1] + \
x[0]*r)*(Power(a,2) + Power(r,2))*(Power(a,2)*Power(x[2],2) + \
Power(r,4))*r_X - 2*r*(Power(a,2) + \
Power(r,2))*(Power(a,2)*Power(x[2],2) + Power(r,4))*(r + x[0]*r_X))) \
+ 4*G*M*Power(r,3)*Power(Power(a,2) + \
Power(r,2),2)*(4*Power(r,4)*(a*x[1] + x[0]*r)*(Power(a,2) + \
Power(r,2))*r_X + 2*Power(r,2)*(a*x[1] + \
x[0]*r)*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_X - 3*(a*x[1] + \
x[0]*r)*(Power(a,2) + Power(r,2))*(Power(a,2)*Power(x[2],2) + \
Power(r,4))*r_X - r*(Power(a,2) + \
Power(r,2))*(Power(a,2)*Power(x[2],2) + Power(r,4))*(r + x[0]*r_X)) + \
2*G*M*x[2]*r*(Power(a,2) + Power(r,2))*(4*Power(r,3)*Power(a*x[1] + \
x[0]*r,2)*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_Z - \
2*x[0]*Power(r,2)*(a*x[1] + x[0]*r)*(Power(a,2) + \
Power(r,2))*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_Z - \
3*r*Power(a*x[1] + x[0]*r,2)*(Power(a,2) + \
Power(r,2))*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_Z + \
2*Power(r,2)*Power(a*x[1] + x[0]*r,2)*(Power(a,2) + \
Power(r,2))*(Power(a,2)*x[2] + 2*Power(r,3)*r_Z) - \
8*x[2]*Power(r,4)*(a*x[1] + x[0]*r)*Power(Power(a,2) + \
Power(r,2),2)*r_X - 4*x[2]*Power(r,2)*(a*x[1] + x[0]*r)*(Power(a,2) + \
Power(r,2))*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_X + \
4*x[2]*(a*x[1] + x[0]*r)*Power(Power(a,2) + \
Power(r,2),2)*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_X + \
2*x[2]*r*Power(Power(a,2) + Power(r,2),2)*(Power(a,2)*Power(x[2],2) + \
Power(r,4))*(r + x[0]*r_X)) - 2*G*M*Power(r,3)*(a*x[0] - \
x[1]*r)*(-2*Power(a,4)*x[1]*Power(x[2],2)*Power(r,3) + \
2*Power(a,2)*x[1]*(-2*Power(a,2) + Power(x[2],2))*Power(r,5) + \
2*x[1]*Power(r,9) - 3*Power(a,6)*x[1]*Power(x[2],2)*(x[1]*r_Y + \
2*x[0]*r_X) - 2*a*Power(r,7)*(a*x[1] - 4*x[0]*x[1]*r_Y - \
4*(Power(x[0],2) - Power(x[1],2))*r_X) - \
4*Power(a,5)*Power(x[2],2)*r*(a*x[1] + 2*x[0]*x[1]*r_Y + \
2*(Power(x[0],2) - Power(x[1],2))*r_X) + \
Power(a,2)*Power(r,4)*(-6*a*x[0]*Power(x[2],2) + \
(Power(a,2)*Power(x[1],2) - Power(x[0],2)*Power(x[2],2))*r_Y + \
2*x[0]*x[1]*(Power(a,2) + Power(x[2],2))*r_X) + \
Power(a,2)*Power(r,6)*(-((Power(x[0],2) - 5*Power(x[1],2))*r_Y) - \
6*x[0]*(a - 2*x[1]*r_X)) + \
Power(a,4)*Power(x[2],2)*Power(r,2)*((-5*Power(x[0],2) + \
Power(x[1],2))*r_Y - 6*x[0]*(a - 2*x[1]*r_X)) + \
3*x[0]*Power(r,8)*(x[0]*r_Y - 2*(a + x[1]*r_X)))) - \
Power(x[5],2)*Power(r,2)*(4*G*M*Power(r,3)*(a*x[1] + \
x[0]*r)*Power(Power(a,2) + Power(r,2),2)*(Power(r,8) - \
3*Power(a,5)*x[0]*Power(x[2],2)*r_Y + \
4*Power(a,4)*x[1]*Power(x[2],2)*r*r_Y + \
2*Power(a,2)*x[1]*Power(x[2],2)*Power(r,3)*r_Y - \
2*x[1]*Power(r,7)*r_Y + Power(a,3)*Power(x[2],2)*Power(r,2)*(a - \
x[0]*r_Y) + a*Power(r,6)*(a + 3*x[0]*r_Y) + \
Power(a,2)*Power(r,4)*(Power(x[2],2) + a*x[0]*r_Y)) - \
2*G*M*Power(r,3)*(a*x[1] + x[0]*r)*Power(a*x[0] - \
x[1]*r,2)*(4*Power(r,4)*(a*x[0] - x[1]*r)*(Power(a,2) + \
Power(r,2))*r_Y + 4*Power(r,2)*(a*x[0] - \
x[1]*r)*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_Y - 3*(a*x[0] - \
x[1]*r)*(Power(a,2) + Power(r,2))*(Power(a,2)*Power(x[2],2) + \
Power(r,4))*r_Y + 2*r*(Power(a,2) + \
Power(r,2))*(Power(a,2)*Power(x[2],2) + Power(r,4))*(r + x[1]*r_Y)) - \
2*G*M*x[2]*r*(a*x[1] + x[0]*r)*(Power(a,2) + \
Power(r,2))*(4*Power(r,3)*Power(a*x[0] - \
x[1]*r,2)*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_Z + \
2*x[1]*Power(r,2)*(a*x[0] - x[1]*r)*(Power(a,2) + \
Power(r,2))*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_Z - \
3*r*Power(a*x[0] - x[1]*r,2)*(Power(a,2) + \
Power(r,2))*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_Z + \
2*Power(r,2)*Power(a*x[0] - x[1]*r,2)*(Power(a,2) + \
Power(r,2))*(Power(a,2)*x[2] + 2*Power(r,3)*r_Z) + \
8*x[2]*Power(r,4)*(a*x[0] - x[1]*r)*Power(Power(a,2) + \
Power(r,2),2)*r_Y + 4*x[2]*Power(r,2)*(a*x[0] - x[1]*r)*(Power(a,2) + \
Power(r,2))*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_Y - \
4*x[2]*(a*x[0] - x[1]*r)*Power(Power(a,2) + \
Power(r,2),2)*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_Y + \
2*x[2]*r*Power(Power(a,2) + Power(r,2),2)*(Power(a,2)*Power(x[2],2) + \
Power(r,4))*(r + x[1]*r_Y)) - (2*G*M*Power(r,3)*Power(a*x[1] + \
x[0]*r,2) - Power(Power(a,2) + \
Power(r,2),2)*(Power(a,2)*Power(x[2],2) + \
Power(r,4)))*(-2*Power(a,4)*x[0]*Power(x[2],2)*Power(r,3) + \
2*Power(a,2)*x[0]*(-2*Power(a,2) + Power(x[2],2))*Power(r,5) + \
2*x[0]*Power(r,9) - 3*Power(a,6)*x[0]*Power(x[2],2)*(2*x[1]*r_Y + \
x[0]*r_X) + 3*x[1]*Power(r,8)*(2*a - 2*x[0]*r_Y + x[1]*r_X) + \
Power(a,4)*Power(x[2],2)*Power(r,2)*(6*a*x[1] + 12*x[0]*x[1]*r_Y + \
(Power(x[0],2) - 5*Power(x[1],2))*r_X) + \
Power(a,2)*Power(r,6)*(6*a*x[1] + 12*x[0]*x[1]*r_Y + (5*Power(x[0],2) \
- Power(x[1],2))*r_X) + Power(a,2)*Power(r,4)*(6*a*x[1]*Power(x[2],2) \
+ 2*x[0]*x[1]*(Power(a,2) + Power(x[2],2))*r_Y + \
(Power(a,2)*Power(x[0],2) - Power(x[1],2)*Power(x[2],2))*r_X) - \
4*Power(a,5)*Power(x[2],2)*r*(2*(Power(x[0],2) - Power(x[1],2))*r_Y + \
x[0]*(a - 2*x[1]*r_X)) - 2*a*Power(r,7)*(-4*(Power(x[0],2) - \
Power(x[1],2))*r_Y + x[0]*(a + 4*x[1]*r_X)))) - \
2*x[4]*x[5]*Power(r,2)*(a*x[1] + \
x[0]*r)*((2*G*M*Power(r,3)*Power(a*x[1] + x[0]*r,2) - \
Power(Power(a,2) + Power(r,2),2)*(Power(a,2)*Power(x[2],2) + \
Power(r,4)))*(4*Power(r,4)*(a*x[1] + x[0]*r)*(Power(a,2) + \
Power(r,2))*r_Y + 4*Power(r,2)*(a*x[1] + \
x[0]*r)*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_Y - 3*(a*x[1] + \
x[0]*r)*(Power(a,2) + Power(r,2))*(Power(a,2)*Power(x[2],2) + \
Power(r,4))*r_Y - 2*r*(Power(a,2) + \
Power(r,2))*(Power(a,2)*Power(x[2],2) + Power(r,4))*(a + x[0]*r_Y)) - \
2*G*M*Power(r,3)*Power(Power(a,2) + \
Power(r,2),2)*((-3*Power(a,5)*x[1]*Power(x[2],2) - \
4*Power(a,4)*x[0]*Power(x[2],2)*r - \
Power(a,3)*x[1]*Power(x[2],2)*Power(r,2) - \
2*Power(a,2)*x[0]*Power(x[2],2)*Power(r,3) + \
Power(a,3)*x[1]*Power(r,4) + 3*a*x[1]*Power(r,6) + \
2*x[0]*Power(r,7))*r_Y + (3*Power(a,5)*x[0]*Power(x[2],2) - \
4*Power(a,4)*x[1]*Power(x[2],2)*r + \
Power(a,3)*x[0]*Power(x[2],2)*Power(r,2) - \
2*Power(a,2)*x[1]*Power(x[2],2)*Power(r,3) - \
Power(a,3)*x[0]*Power(r,4) - 3*a*x[0]*Power(r,6) + \
2*x[1]*Power(r,7))*r_X) - 2*G*M*Power(r,3)*Power(a*x[0] - \
x[1]*r,2)*(4*Power(r,4)*(a*x[0] - x[1]*r)*(Power(a,2) + \
Power(r,2))*r_X + 4*Power(r,2)*(a*x[0] - \
x[1]*r)*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_X - 3*(a*x[0] - \
x[1]*r)*(Power(a,2) + Power(r,2))*(Power(a,2)*Power(x[2],2) + \
Power(r,4))*r_X - 2*r*(Power(a,2) + \
Power(r,2))*(Power(a,2)*Power(x[2],2) + Power(r,4))*(a - x[1]*r_X)) + \
2*G*M*x[2]*r*(Power(a,2) + \
Power(r,2))*(2*Power(a,5)*x[2]*Power(r,4)*(x[1]*r_Y - x[0]*r_X) + \
2*Power(a,7)*Power(x[2],3)*(-(x[1]*r_Y) + x[0]*r_X) - \
3*Power(a,6)*Power(x[2],2)*r*(x[0]*x[1]*r_Z + x[0]*x[2]*r_Y + \
x[1]*x[2]*r_X) + 2*Power(a,5)*x[2]*Power(r,2)*(a*x[0]*x[1] + \
2*(-Power(x[0],2) + Power(x[1],2))*x[2]*r_Z - x[1]*Power(x[2],2)*r_Y \
+ x[0]*Power(x[2],2)*r_X) + \
2*Power(a,4)*x[2]*Power(r,3)*(a*Power(x[0],2) - a*Power(x[1],2) + \
3*x[0]*x[1]*x[2]*r_Z - 2*x[0]*Power(x[2],2)*r_Y - \
2*x[1]*Power(x[2],2)*r_X) - \
2*Power(a,2)*x[2]*Power(r,6)*(-3*a*x[1]*r_Y + x[0]*(x[1] + 3*a*r_X)) \
+ 4*a*Power(r,8)*((Power(x[0],2) - Power(x[1],2))*r_Z + \
x[2]*(x[1]*r_Y - x[0]*r_X)) + \
2*Power(a,2)*Power(r,7)*(3*x[0]*x[1]*r_Z + 2*x[2]*(x[0]*r_Y + \
x[1]*r_X)) + Power(r,9)*(-3*x[0]*x[1]*r_Z + 3*x[2]*(x[0]*r_Y + \
x[1]*r_X)) + Power(a,2)*Power(r,5)*(x[0]*x[1]*(Power(a,2) + \
Power(x[2],2))*r_Z + x[2]*(2*a*(Power(x[0],2) - Power(x[1],2)) + \
x[0]*(Power(a,2) - Power(x[2],2))*r_Y + x[1]*(Power(a,2) - \
Power(x[2],2))*r_X)))) + 2*x[5]*x[6]*r*(2*G*M*Power(r,4)*(a*x[1] + \
x[0]*r)*Power(a*x[0] - x[1]*r,2)*(4*Power(r,2)*(a*x[0] - \
x[1]*r)*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_Z + \
2*x[1]*r*(Power(a,2) + Power(r,2))*(Power(a,2)*Power(x[2],2) + \
Power(r,4))*r_Z - 3*(a*x[0] - x[1]*r)*(Power(a,2) + \
Power(r,2))*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_Z + r*(a*x[0] - \
x[1]*r)*(Power(a,2) + Power(r,2))*(2*Power(a,2)*x[2] + \
4*Power(r,3)*r_Z)) + 2*G*M*Power(x[2],3)*r*(a*x[1] + \
x[0]*r)*Power(Power(a,2) + Power(r,2),4)*(Power(a,2)*Power(x[2],2) - \
3*Power(r,4))*r_Y - 2*G*M*Power(r,3)*(a*x[1] + \
x[0]*r)*Power(Power(a,2) + Power(r,2),2)*(2*Power(r,3)*(a*x[0] - \
x[1]*r)*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_Z + \
x[1]*Power(r,2)*(Power(a,2) + Power(r,2))*(Power(a,2)*Power(x[2],2) + \
Power(r,4))*r_Z + 3*r*(-(a*x[0]) + x[1]*r)*(Power(a,2) + \
Power(r,2))*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_Z + \
Power(r,2)*(a*x[0] - x[1]*r)*(Power(a,2) + \
Power(r,2))*(2*Power(a,2)*x[2] + 4*Power(r,3)*r_Z) - \
4*x[2]*Power(r,4)*Power(Power(a,2) + Power(r,2),2)*r_Y + \
2*x[2]*Power(Power(a,2) + Power(r,2),2)*(Power(a,2)*Power(x[2],2) + \
Power(r,4))*r_Y) + (2*G*M*Power(r,3)*Power(a*x[1] + x[0]*r,2) - \
Power(Power(a,2) + Power(r,2),2)*(Power(a,2)*Power(x[2],2) + \
Power(r,4)))*(2*Power(a,7)*Power(x[2],3)*(x[1]*r_Y + x[0]*r_X) - \
2*Power(a,5)*x[2]*Power(r,4)*(x[1]*r_Y + x[0]*r_X) + \
2*Power(a,5)*x[2]*Power(r,2)*(a*x[0]*x[1] + 2*(-Power(x[0],2) + \
Power(x[1],2))*x[2]*r_Z + x[1]*Power(x[2],2)*r_Y + \
x[0]*Power(x[2],2)*r_X) + \
2*Power(a,4)*x[2]*Power(r,3)*(a*Power(x[0],2) - a*Power(x[1],2) + \
2*a*Power(x[2],2) + 3*x[0]*x[1]*x[2]*r_Z + 2*x[0]*Power(x[2],2)*r_Y - \
2*x[1]*Power(x[2],2)*r_X) - \
2*Power(a,2)*x[2]*Power(r,6)*(3*a*x[1]*r_Y + x[0]*(x[1] + 3*a*r_X)) + \
4*a*Power(r,8)*((Power(x[0],2) - Power(x[1],2))*r_Z - x[2]*(x[1]*r_Y \
+ x[0]*r_X)) + Power(a,6)*Power(x[2],2)*r*(-3*x[0]*x[1]*r_Z + \
x[2]*(2*a + 3*x[0]*r_Y - 3*x[1]*r_X)) + \
2*Power(a,2)*Power(r,7)*(3*x[0]*x[1]*r_Z + 2*x[2]*(a - x[0]*r_Y + \
x[1]*r_X)) + Power(r,9)*(-3*x[0]*x[1]*r_Z + x[2]*(2*a - 3*x[0]*r_Y + \
3*x[1]*r_X)) + Power(a,2)*Power(r,5)*(x[0]*x[1]*(Power(a,2) + \
Power(x[2],2))*r_Z + x[2]*(2*a*(Power(a,2) + Power(x[0],2) - \
Power(x[1],2) + Power(x[2],2)) + x[0]*(-Power(a,2) + \
Power(x[2],2))*r_Y + x[1]*(Power(a,2) - Power(x[2],2))*r_X)))) - \
2*x[4]*x[6]*Power(r,2)*(a*x[1] + \
x[0]*r)*((2*G*M*Power(r,3)*Power(a*x[1] + x[0]*r,2) - \
Power(Power(a,2) + Power(r,2),2)*(Power(a,2)*Power(x[2],2) + \
Power(r,4)))*(4*Power(r,2)*(a*x[1] + \
x[0]*r)*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_Z - \
2*x[0]*r*(Power(a,2) + Power(r,2))*(Power(a,2)*Power(x[2],2) + \
Power(r,4))*r_Z - 3*(a*x[1] + x[0]*r)*(Power(a,2) + \
Power(r,2))*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_Z + 2*r*(a*x[1] \
+ x[0]*r)*(Power(a,2) + Power(r,2))*(Power(a,2)*x[2] + \
2*Power(r,3)*r_Z)) - 2*G*M*Power(x[2],3)*Power(Power(a,2) + \
Power(r,2),4)*(Power(a,2)*Power(x[2],2) - 3*Power(r,4))*r_X - \
2*G*M*Power(r,2)*Power(Power(a,2) + \
Power(r,2),2)*(2*Power(r,3)*(a*x[1] + \
x[0]*r)*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_Z - \
x[0]*Power(r,2)*(Power(a,2) + Power(r,2))*(Power(a,2)*Power(x[2],2) + \
Power(r,4))*r_Z - 3*r*(a*x[1] + x[0]*r)*(Power(a,2) + \
Power(r,2))*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_Z + \
2*Power(r,2)*(a*x[1] + x[0]*r)*(Power(a,2) + \
Power(r,2))*(Power(a,2)*x[2] + 2*Power(r,3)*r_Z) + \
4*x[2]*Power(r,4)*Power(Power(a,2) + Power(r,2),2)*r_X - \
2*x[2]*Power(Power(a,2) + Power(r,2),2)*(Power(a,2)*Power(x[2],2) + \
Power(r,4))*r_X) + 2*G*M*Power(r,2)*(a*x[0] - \
x[1]*r)*(-2*Power(a,7)*Power(x[2],3)*(x[1]*r_Y + x[0]*r_X) + \
2*Power(a,5)*x[2]*Power(r,4)*(x[1]*r_Y + x[0]*r_X) + \
2*Power(a,2)*x[2]*Power(r,6)*(-(x[0]*x[1]) + 3*a*x[1]*r_Y + \
3*a*x[0]*r_X) + 2*Power(a,5)*x[2]*Power(r,2)*(a*x[0]*x[1] + \
2*(-Power(x[0],2) + Power(x[1],2))*x[2]*r_Z - x[1]*Power(x[2],2)*r_Y \
- x[0]*Power(x[2],2)*r_X) + \
2*Power(a,4)*x[2]*Power(r,3)*(a*Power(x[0],2) - a*Power(x[1],2) - \
2*a*Power(x[2],2) + 3*x[0]*x[1]*x[2]*r_Z - 2*x[0]*Power(x[2],2)*r_Y + \
2*x[1]*Power(x[2],2)*r_X) + 4*a*Power(r,8)*((Power(x[0],2) - \
Power(x[1],2))*r_Z + x[2]*(x[1]*r_Y + x[0]*r_X)) - \
Power(a,6)*Power(x[2],2)*r*(3*x[0]*x[1]*r_Z + x[2]*(2*a + 3*x[0]*r_Y \
- 3*x[1]*r_X)) + 2*Power(a,2)*Power(r,7)*(3*x[0]*x[1]*r_Z - 2*x[2]*(a \
- x[0]*r_Y + x[1]*r_X)) - Power(r,9)*(3*x[0]*x[1]*r_Z + x[2]*(2*a - \
3*x[0]*r_Y + 3*x[1]*r_X)) + \
Power(a,2)*Power(r,5)*(x[0]*x[1]*(Power(a,2) + Power(x[2],2))*r_Z - \
x[2]*(2*a*(Power(a,2) - Power(x[0],2) + Power(x[1],2) + \
Power(x[2],2)) + x[0]*(-Power(a,2) + Power(x[2],2))*r_Y + \
x[1]*(Power(a,2) - Power(x[2],2))*r_X))))))/(Power(Power(a,2) + \
Power(r,2),5)*Power(Power(a,2)*Power(x[2],2) + Power(r,4),3));
        dxdt[5] = (G*M*(Power(x[7],2)*Power(Power(a,2) + \
Power(r,2),3)*(2*G*M*x[2]*Power(r,4)*(-(a*x[0]) + x[1]*r)*(Power(a,2) \
+ Power(r,2))*(2*Power(a,2)*x[2]*r - 3*Power(a,2)*Power(x[2],2)*r_Z + \
Power(r,4)*r_Z) + (-3*Power(a,2)*Power(x[2],2)*Power(r,2) + \
Power(r,6))*(2*G*M*Power(r,3)*Power(a*x[0] - x[1]*r,2) - \
Power(Power(a,2) + Power(r,2),2)*(Power(a,2)*Power(x[2],2) + \
Power(r,4)))*r_Y + 2*G*M*Power(r,3)*(a*x[1] + x[0]*r)*(-(a*x[0]) + \
x[1]*r)*(-3*Power(a,2)*Power(x[2],2)*Power(r,2) + Power(r,6))*r_X) + \
2*x[7]*x[6]*r*(Power(a,2) + Power(r,2))*(2*G*M*Power(r,4)*(-(a*x[0]) \
+ x[1]*r)*Power(Power(a,2) + Power(r,2),3)*(2*Power(a,2)*x[2]*r - \
3*Power(a,2)*Power(x[2],2)*r_Z + Power(r,4)*r_Z) + \
(2*G*M*Power(r,3)*Power(a*x[0] - x[1]*r,2) - Power(Power(a,2) + \
Power(r,2),2)*(Power(a,2)*Power(x[2],2) + \
Power(r,4)))*(2*Power(r,3)*(a*x[0] - \
x[1]*r)*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_Z + \
x[1]*Power(r,2)*(Power(a,2) + Power(r,2))*(Power(a,2)*Power(x[2],2) + \
Power(r,4))*r_Z + 3*r*(-(a*x[0]) + x[1]*r)*(Power(a,2) + \
Power(r,2))*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_Z + \
Power(r,2)*(a*x[0] - x[1]*r)*(Power(a,2) + \
Power(r,2))*(2*Power(a,2)*x[2] + 4*Power(r,3)*r_Z) + \
4*x[2]*Power(r,4)*Power(Power(a,2) + Power(r,2),2)*r_Y - \
2*x[2]*Power(Power(a,2) + Power(r,2),2)*(Power(a,2)*Power(x[2],2) + \
Power(r,4))*r_Y) + 2*G*M*Power(r,3)*(a*x[1] + x[0]*r)*(a*x[0] - \
x[1]*r)*(2*Power(r,3)*(a*x[1] + x[0]*r)*(Power(a,2)*Power(x[2],2) + \
Power(r,4))*r_Z - x[0]*Power(r,2)*(Power(a,2) + \
Power(r,2))*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_Z - 3*r*(a*x[1] \
+ x[0]*r)*(Power(a,2) + Power(r,2))*(Power(a,2)*Power(x[2],2) + \
Power(r,4))*r_Z + 2*Power(r,2)*(a*x[1] + x[0]*r)*(Power(a,2) + \
Power(r,2))*(Power(a,2)*x[2] + 2*Power(r,3)*r_Z) - \
4*x[2]*Power(r,4)*Power(Power(a,2) + Power(r,2),2)*r_X + \
2*x[2]*Power(Power(a,2) + Power(r,2),2)*(Power(a,2)*Power(x[2],2) + \
Power(r,4))*r_X)) - Power(x[6],2)*(Power(a,2) + \
Power(r,2))*(4*G*M*Power(r,4)*(-(a*x[0]) + x[1]*r)*Power(Power(a,2) + \
Power(r,2),3)*(-(Power(a,2)*Power(x[2],2)) + Power(r,4))*(r - \
2*x[2]*r_Z) - 2*G*M*Power(x[2],2)*Power(r,2)*(-(a*x[0]) + \
x[1]*r)*Power(Power(a,2) + Power(r,2),3)*(2*Power(r,5) + \
Power(a,2)*Power(x[2],3)*r_Z - 3*x[2]*Power(r,4)*r_Z) + \
(2*G*M*Power(r,3)*Power(a*x[0] - x[1]*r,2) - Power(Power(a,2) + \
Power(r,2),2)*(Power(a,2)*Power(x[2],2) + \
Power(r,4)))*(2*Power(a,4)*x[1]*Power(x[2],2)*Power(r,3) - \
2*x[1]*Power(r,9) + 4*Power(a,5)*x[0]*Power(x[2],3)*r*r_Z + \
2*Power(a,2)*x[2]*Power(r,5)*(x[1]*x[2] - 2*a*x[0]*r_Z) - \
2*a*Power(r,7)*(a*x[1] + 4*x[0]*x[2]*r_Z) + \
Power(a,6)*Power(x[2],4)*r_Y + 2*Power(a,2)*Power(r,6)*(a*x[0] + \
x[1]*x[2]*r_Z - 3*Power(x[2],2)*r_Y) + Power(r,8)*(2*a*x[0] + \
6*x[1]*x[2]*r_Z - 3*Power(x[2],2)*r_Y) - \
2*Power(a,4)*Power(x[2],2)*Power(r,2)*(a*x[0] + 3*x[1]*x[2]*r_Z - \
Power(x[2],2)*r_Y) + Power(a,2)*Power(x[2],2)*Power(r,4)*(-2*a*x[0] - \
2*x[1]*x[2]*r_Z + (-3*Power(a,2) + Power(x[2],2))*r_Y)) + \
2*G*M*Power(r,3)*(a*x[1] + x[0]*r)*(a*x[0] - \
x[1]*r)*(2*Power(r,2)*(a*x[1] + x[0]*r)*(Power(a,2) + \
Power(r,2))*(Power(a,2)*Power(x[2],2) + Power(r,4)) - \
4*x[2]*Power(r,3)*(a*x[1] + x[0]*r)*(Power(a,2)*Power(x[2],2) + \
Power(r,4))*r_Z + 2*x[0]*x[2]*Power(r,2)*(Power(a,2) + \
Power(r,2))*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_Z + \
4*x[2]*r*(a*x[1] + x[0]*r)*(Power(a,2) + \
Power(r,2))*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_Z - \
4*x[2]*Power(r,2)*(a*x[1] + x[0]*r)*(Power(a,2) + \
Power(r,2))*(Power(a,2)*x[2] + 2*Power(r,3)*r_Z) + \
4*Power(x[2],2)*Power(r,4)*Power(Power(a,2) + Power(r,2),2)*r_X - \
Power(x[2],2)*Power(Power(a,2) + \
Power(r,2),2)*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_X)) - \
4*G*M*x[7]*x[5]*Power(r,3)*(a*x[0] - x[1]*r)*(Power(a,2) + \
Power(r,2))*(Power(Power(a,2) + \
Power(r,2),3)*(-3*Power(a,2)*Power(x[2],2)*Power(r,2) + \
Power(r,6))*r_Y - x[2]*(Power(a,2) + \
Power(r,2))*(2*Power(r,3)*(a*x[0] - x[1]*r)*(Power(a,2)*Power(x[2],2) \
+ Power(r,4))*r_Z + x[1]*Power(r,2)*(Power(a,2) + \
Power(r,2))*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_Z + \
3*r*(-(a*x[0]) + x[1]*r)*(Power(a,2) + \
Power(r,2))*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_Z + \
Power(r,2)*(a*x[0] - x[1]*r)*(Power(a,2) + \
Power(r,2))*(2*Power(a,2)*x[2] + 4*Power(r,3)*r_Z) + \
4*x[2]*Power(r,4)*Power(Power(a,2) + Power(r,2),2)*r_Y - \
2*x[2]*Power(Power(a,2) + Power(r,2),2)*(Power(a,2)*Power(x[2],2) + \
Power(r,4))*r_Y) + Power(r,2)*(a*x[1] + \
x[0]*r)*(2*Power(a,3)*Power(r,5) + \
3*Power(a,5)*Power(x[2],2)*(x[1]*r_Y + x[0]*r_X) + \
Power(a,3)*Power(x[2],2)*Power(r,2)*(x[1]*r_Y + x[0]*r_X) - \
Power(a,3)*Power(r,4)*(x[1]*r_Y + x[0]*r_X) - \
3*a*Power(r,6)*(x[1]*r_Y + x[0]*r_X) + \
2*Power(a,4)*Power(x[2],2)*r*(a + 2*x[0]*r_Y - 2*x[1]*r_X) + \
2*Power(a,2)*Power(x[2],2)*Power(r,3)*(a + x[0]*r_Y - x[1]*r_X) + \
2*Power(r,7)*(a - x[0]*r_Y + x[1]*r_X))) - \
2*x[7]*x[4]*Power(r,2)*(Power(a,2) + Power(r,2))*(-2*G*M*r*(-(a*x[0]) \
+ x[1]*r)*Power(Power(a,2) + \
Power(r,2),3)*(-3*Power(a,2)*Power(x[2],2)*Power(r,2) + \
Power(r,6))*r_X + 2*G*M*x[2]*r*(a*x[0] - x[1]*r)*(Power(a,2) + \
Power(r,2))*(2*Power(r,3)*(a*x[1] + x[0]*r)*(Power(a,2)*Power(x[2],2) \
+ Power(r,4))*r_Z - x[0]*Power(r,2)*(Power(a,2) + \
Power(r,2))*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_Z - 3*r*(a*x[1] \
+ x[0]*r)*(Power(a,2) + Power(r,2))*(Power(a,2)*Power(x[2],2) + \
Power(r,4))*r_Z + 2*Power(r,2)*(a*x[1] + x[0]*r)*(Power(a,2) + \
Power(r,2))*(Power(a,2)*x[2] + 2*Power(r,3)*r_Z) - \
4*x[2]*Power(r,4)*Power(Power(a,2) + Power(r,2),2)*r_X + \
2*x[2]*Power(Power(a,2) + Power(r,2),2)*(Power(a,2)*Power(x[2],2) + \
Power(r,4))*r_X) + (2*G*M*Power(r,3)*Power(a*x[0] - x[1]*r,2) - \
Power(Power(a,2) + Power(r,2),2)*(Power(a,2)*Power(x[2],2) + \
Power(r,4)))*(2*Power(a,3)*Power(r,5) + \
3*Power(a,5)*Power(x[2],2)*(x[1]*r_Y + x[0]*r_X) + \
Power(a,3)*Power(x[2],2)*Power(r,2)*(x[1]*r_Y + x[0]*r_X) - \
Power(a,3)*Power(r,4)*(x[1]*r_Y + x[0]*r_X) - \
3*a*Power(r,6)*(x[1]*r_Y + x[0]*r_X) + \
2*Power(a,4)*Power(x[2],2)*r*(a + 2*x[0]*r_Y - 2*x[1]*r_X) + \
2*Power(a,2)*Power(x[2],2)*Power(r,3)*(a + x[0]*r_Y - x[1]*r_X) + \
2*Power(r,7)*(a - x[0]*r_Y + x[1]*r_X))) - \
Power(x[4],2)*Power(r,2)*(-2*G*M*Power(r,3)*Power(a*x[1] + \
x[0]*r,2)*(a*x[0] - x[1]*r)*(4*Power(r,4)*(a*x[1] + \
x[0]*r)*(Power(a,2) + Power(r,2))*r_X + 4*Power(r,2)*(a*x[1] + \
x[0]*r)*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_X - 3*(a*x[1] + \
x[0]*r)*(Power(a,2) + Power(r,2))*(Power(a,2)*Power(x[2],2) + \
Power(r,4))*r_X - 2*r*(Power(a,2) + \
Power(r,2))*(Power(a,2)*Power(x[2],2) + Power(r,4))*(r + x[0]*r_X)) + \
4*G*M*Power(r,3)*(a*x[0] - x[1]*r)*Power(Power(a,2) + \
Power(r,2),2)*(4*Power(r,4)*(a*x[1] + x[0]*r)*(Power(a,2) + \
Power(r,2))*r_X + 2*Power(r,2)*(a*x[1] + \
x[0]*r)*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_X - 3*(a*x[1] + \
x[0]*r)*(Power(a,2) + Power(r,2))*(Power(a,2)*Power(x[2],2) + \
Power(r,4))*r_X - r*(Power(a,2) + \
Power(r,2))*(Power(a,2)*Power(x[2],2) + Power(r,4))*(r + x[0]*r_X)) + \
2*G*M*x[2]*r*(a*x[0] - x[1]*r)*(Power(a,2) + \
Power(r,2))*(4*Power(r,3)*Power(a*x[1] + \
x[0]*r,2)*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_Z - \
2*x[0]*Power(r,2)*(a*x[1] + x[0]*r)*(Power(a,2) + \
Power(r,2))*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_Z - \
3*r*Power(a*x[1] + x[0]*r,2)*(Power(a,2) + \
Power(r,2))*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_Z + \
2*Power(r,2)*Power(a*x[1] + x[0]*r,2)*(Power(a,2) + \
Power(r,2))*(Power(a,2)*x[2] + 2*Power(r,3)*r_Z) - \
8*x[2]*Power(r,4)*(a*x[1] + x[0]*r)*Power(Power(a,2) + \
Power(r,2),2)*r_X - 4*x[2]*Power(r,2)*(a*x[1] + x[0]*r)*(Power(a,2) + \
Power(r,2))*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_X + \
4*x[2]*(a*x[1] + x[0]*r)*Power(Power(a,2) + \
Power(r,2),2)*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_X + \
2*x[2]*r*Power(Power(a,2) + Power(r,2),2)*(Power(a,2)*Power(x[2],2) + \
Power(r,4))*(r + x[0]*r_X)) - (2*G*M*Power(r,3)*Power(a*x[0] - \
x[1]*r,2) - Power(Power(a,2) + \
Power(r,2),2)*(Power(a,2)*Power(x[2],2) + \
Power(r,4)))*(-2*Power(a,4)*x[1]*Power(x[2],2)*Power(r,3) + \
2*Power(a,2)*x[1]*(-2*Power(a,2) + Power(x[2],2))*Power(r,5) + \
2*x[1]*Power(r,9) - 3*Power(a,6)*x[1]*Power(x[2],2)*(x[1]*r_Y + \
2*x[0]*r_X) - 2*a*Power(r,7)*(a*x[1] - 4*x[0]*x[1]*r_Y - \
4*(Power(x[0],2) - Power(x[1],2))*r_X) - \
4*Power(a,5)*Power(x[2],2)*r*(a*x[1] + 2*x[0]*x[1]*r_Y + \
2*(Power(x[0],2) - Power(x[1],2))*r_X) + \
Power(a,2)*Power(r,4)*(-6*a*x[0]*Power(x[2],2) + \
(Power(a,2)*Power(x[1],2) - Power(x[0],2)*Power(x[2],2))*r_Y + \
2*x[0]*x[1]*(Power(a,2) + Power(x[2],2))*r_X) + \
Power(a,2)*Power(r,6)*(-((Power(x[0],2) - 5*Power(x[1],2))*r_Y) - \
6*x[0]*(a - 2*x[1]*r_X)) + \
Power(a,4)*Power(x[2],2)*Power(r,2)*((-5*Power(x[0],2) + \
Power(x[1],2))*r_Y - 6*x[0]*(a - 2*x[1]*r_X)) + \
3*x[0]*Power(r,8)*(x[0]*r_Y - 2*(a + x[1]*r_X)))) + \
Power(x[5],2)*Power(r,2)*(a*x[0] - \
x[1]*r)*(4*G*M*Power(r,3)*Power(Power(a,2) + \
Power(r,2),2)*(Power(r,8) - 3*Power(a,5)*x[0]*Power(x[2],2)*r_Y + \
4*Power(a,4)*x[1]*Power(x[2],2)*r*r_Y + \
2*Power(a,2)*x[1]*Power(x[2],2)*Power(r,3)*r_Y - \
2*x[1]*Power(r,7)*r_Y + Power(a,3)*Power(x[2],2)*Power(r,2)*(a - \
x[0]*r_Y) + a*Power(r,6)*(a + 3*x[0]*r_Y) + \
Power(a,2)*Power(r,4)*(Power(x[2],2) + a*x[0]*r_Y)) - \
(2*G*M*Power(r,3)*Power(a*x[0] - x[1]*r,2) - Power(Power(a,2) + \
Power(r,2),2)*(Power(a,2)*Power(x[2],2) + \
Power(r,4)))*(4*Power(r,4)*(a*x[0] - x[1]*r)*(Power(a,2) + \
Power(r,2))*r_Y + 4*Power(r,2)*(a*x[0] - \
x[1]*r)*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_Y - 3*(a*x[0] - \
x[1]*r)*(Power(a,2) + Power(r,2))*(Power(a,2)*Power(x[2],2) + \
Power(r,4))*r_Y + 2*r*(Power(a,2) + \
Power(r,2))*(Power(a,2)*Power(x[2],2) + Power(r,4))*(r + x[1]*r_Y)) - \
2*G*M*x[2]*r*(Power(a,2) + Power(r,2))*(4*Power(r,3)*Power(a*x[0] - \
x[1]*r,2)*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_Z + \
2*x[1]*Power(r,2)*(a*x[0] - x[1]*r)*(Power(a,2) + \
Power(r,2))*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_Z - \
3*r*Power(a*x[0] - x[1]*r,2)*(Power(a,2) + \
Power(r,2))*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_Z + \
2*Power(r,2)*Power(a*x[0] - x[1]*r,2)*(Power(a,2) + \
Power(r,2))*(Power(a,2)*x[2] + 2*Power(r,3)*r_Z) + \
8*x[2]*Power(r,4)*(a*x[0] - x[1]*r)*Power(Power(a,2) + \
Power(r,2),2)*r_Y + 4*x[2]*Power(r,2)*(a*x[0] - x[1]*r)*(Power(a,2) + \
Power(r,2))*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_Y - \
4*x[2]*(a*x[0] - x[1]*r)*Power(Power(a,2) + \
Power(r,2),2)*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_Y + \
2*x[2]*r*Power(Power(a,2) + Power(r,2),2)*(Power(a,2)*Power(x[2],2) + \
Power(r,4))*(r + x[1]*r_Y)) - 2*G*M*Power(r,3)*(a*x[1] + \
x[0]*r)*(-2*Power(a,4)*x[0]*Power(x[2],2)*Power(r,3) + \
2*Power(a,2)*x[0]*(-2*Power(a,2) + Power(x[2],2))*Power(r,5) + \
2*x[0]*Power(r,9) - 3*Power(a,6)*x[0]*Power(x[2],2)*(2*x[1]*r_Y + \
x[0]*r_X) + 3*x[1]*Power(r,8)*(2*a - 2*x[0]*r_Y + x[1]*r_X) + \
Power(a,4)*Power(x[2],2)*Power(r,2)*(6*a*x[1] + 12*x[0]*x[1]*r_Y + \
(Power(x[0],2) - 5*Power(x[1],2))*r_X) + \
Power(a,2)*Power(r,6)*(6*a*x[1] + 12*x[0]*x[1]*r_Y + (5*Power(x[0],2) \
- Power(x[1],2))*r_X) + Power(a,2)*Power(r,4)*(6*a*x[1]*Power(x[2],2) \
+ 2*x[0]*x[1]*(Power(a,2) + Power(x[2],2))*r_Y + \
(Power(a,2)*Power(x[0],2) - Power(x[1],2)*Power(x[2],2))*r_X) - \
4*Power(a,5)*Power(x[2],2)*r*(2*(Power(x[0],2) - Power(x[1],2))*r_Y + \
x[0]*(a - 2*x[1]*r_X)) - 2*a*Power(r,7)*(-4*(Power(x[0],2) - \
Power(x[1],2))*r_Y + x[0]*(a + 4*x[1]*r_X)))) + \
2*x[4]*x[5]*Power(r,2)*(2*G*M*Power(r,3)*Power(a*x[1] + \
x[0]*r,2)*(a*x[0] - x[1]*r)*(4*Power(r,4)*(a*x[1] + \
x[0]*r)*(Power(a,2) + Power(r,2))*r_Y + 4*Power(r,2)*(a*x[1] + \
x[0]*r)*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_Y - 3*(a*x[1] + \
x[0]*r)*(Power(a,2) + Power(r,2))*(Power(a,2)*Power(x[2],2) + \
Power(r,4))*r_Y - 2*r*(Power(a,2) + \
Power(r,2))*(Power(a,2)*Power(x[2],2) + Power(r,4))*(a + x[0]*r_Y)) + \
2*G*M*Power(r,3)*(-(a*x[0]) + x[1]*r)*Power(Power(a,2) + \
Power(r,2),2)*((-3*Power(a,5)*x[1]*Power(x[2],2) - \
4*Power(a,4)*x[0]*Power(x[2],2)*r - \
Power(a,3)*x[1]*Power(x[2],2)*Power(r,2) - \
2*Power(a,2)*x[0]*Power(x[2],2)*Power(r,3) + \
Power(a,3)*x[1]*Power(r,4) + 3*a*x[1]*Power(r,6) + \
2*x[0]*Power(r,7))*r_Y + (3*Power(a,5)*x[0]*Power(x[2],2) - \
4*Power(a,4)*x[1]*Power(x[2],2)*r + \
Power(a,3)*x[0]*Power(x[2],2)*Power(r,2) - \
2*Power(a,2)*x[1]*Power(x[2],2)*Power(r,3) - \
Power(a,3)*x[0]*Power(r,4) - 3*a*x[0]*Power(r,6) + \
2*x[1]*Power(r,7))*r_X) - (a*x[0] - \
x[1]*r)*(2*G*M*Power(r,3)*Power(a*x[0] - x[1]*r,2) - Power(Power(a,2) \
+ Power(r,2),2)*(Power(a,2)*Power(x[2],2) + \
Power(r,4)))*(4*Power(r,4)*(a*x[0] - x[1]*r)*(Power(a,2) + \
Power(r,2))*r_X + 4*Power(r,2)*(a*x[0] - \
x[1]*r)*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_X - 3*(a*x[0] - \
x[1]*r)*(Power(a,2) + Power(r,2))*(Power(a,2)*Power(x[2],2) + \
Power(r,4))*r_X - 2*r*(Power(a,2) + \
Power(r,2))*(Power(a,2)*Power(x[2],2) + Power(r,4))*(a - x[1]*r_X)) + \
2*G*M*x[2]*r*(a*x[0] - x[1]*r)*(Power(a,2) + \
Power(r,2))*(2*Power(a,5)*x[2]*Power(r,4)*(x[1]*r_Y - x[0]*r_X) + \
2*Power(a,7)*Power(x[2],3)*(-(x[1]*r_Y) + x[0]*r_X) - \
3*Power(a,6)*Power(x[2],2)*r*(x[0]*x[1]*r_Z + x[0]*x[2]*r_Y + \
x[1]*x[2]*r_X) + 2*Power(a,5)*x[2]*Power(r,2)*(a*x[0]*x[1] + \
2*(-Power(x[0],2) + Power(x[1],2))*x[2]*r_Z - x[1]*Power(x[2],2)*r_Y \
+ x[0]*Power(x[2],2)*r_X) + \
2*Power(a,4)*x[2]*Power(r,3)*(a*Power(x[0],2) - a*Power(x[1],2) + \
3*x[0]*x[1]*x[2]*r_Z - 2*x[0]*Power(x[2],2)*r_Y - \
2*x[1]*Power(x[2],2)*r_X) - \
2*Power(a,2)*x[2]*Power(r,6)*(-3*a*x[1]*r_Y + x[0]*(x[1] + 3*a*r_X)) \
+ 4*a*Power(r,8)*((Power(x[0],2) - Power(x[1],2))*r_Z + \
x[2]*(x[1]*r_Y - x[0]*r_X)) + \
2*Power(a,2)*Power(r,7)*(3*x[0]*x[1]*r_Z + 2*x[2]*(x[0]*r_Y + \
x[1]*r_X)) + Power(r,9)*(-3*x[0]*x[1]*r_Z + 3*x[2]*(x[0]*r_Y + \
x[1]*r_X)) + Power(a,2)*Power(r,5)*(x[0]*x[1]*(Power(a,2) + \
Power(x[2],2))*r_Z + x[2]*(2*a*(Power(x[0],2) - Power(x[1],2)) + \
x[0]*(Power(a,2) - Power(x[2],2))*r_Y + x[1]*(Power(a,2) - \
Power(x[2],2))*r_X)))) - 2*x[5]*x[6]*Power(r,2)*(a*x[0] - \
x[1]*r)*((2*G*M*Power(r,3)*Power(a*x[0] - x[1]*r,2) - \
Power(Power(a,2) + Power(r,2),2)*(Power(a,2)*Power(x[2],2) + \
Power(r,4)))*(4*Power(r,2)*(a*x[0] - \
x[1]*r)*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_Z + \
2*x[1]*r*(Power(a,2) + Power(r,2))*(Power(a,2)*Power(x[2],2) + \
Power(r,4))*r_Z - 3*(a*x[0] - x[1]*r)*(Power(a,2) + \
Power(r,2))*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_Z + r*(a*x[0] - \
x[1]*r)*(Power(a,2) + Power(r,2))*(2*Power(a,2)*x[2] + \
4*Power(r,3)*r_Z)) + 2*G*M*Power(x[2],3)*Power(Power(a,2) + \
Power(r,2),4)*(Power(a,2)*Power(x[2],2) - 3*Power(r,4))*r_Y - \
2*G*M*Power(r,2)*Power(Power(a,2) + \
Power(r,2),2)*(2*Power(r,3)*(a*x[0] - \
x[1]*r)*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_Z + \
x[1]*Power(r,2)*(Power(a,2) + Power(r,2))*(Power(a,2)*Power(x[2],2) + \
Power(r,4))*r_Z + 3*r*(-(a*x[0]) + x[1]*r)*(Power(a,2) + \
Power(r,2))*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_Z + \
Power(r,2)*(a*x[0] - x[1]*r)*(Power(a,2) + \
Power(r,2))*(2*Power(a,2)*x[2] + 4*Power(r,3)*r_Z) - \
4*x[2]*Power(r,4)*Power(Power(a,2) + Power(r,2),2)*r_Y + \
2*x[2]*Power(Power(a,2) + Power(r,2),2)*(Power(a,2)*Power(x[2],2) + \
Power(r,4))*r_Y) + 2*G*M*Power(r,2)*(a*x[1] + \
x[0]*r)*(2*Power(a,7)*Power(x[2],3)*(x[1]*r_Y + x[0]*r_X) - \
2*Power(a,5)*x[2]*Power(r,4)*(x[1]*r_Y + x[0]*r_X) + \
2*Power(a,5)*x[2]*Power(r,2)*(a*x[0]*x[1] + 2*(-Power(x[0],2) + \
Power(x[1],2))*x[2]*r_Z + x[1]*Power(x[2],2)*r_Y + \
x[0]*Power(x[2],2)*r_X) + \
2*Power(a,4)*x[2]*Power(r,3)*(a*Power(x[0],2) - a*Power(x[1],2) + \
2*a*Power(x[2],2) + 3*x[0]*x[1]*x[2]*r_Z + 2*x[0]*Power(x[2],2)*r_Y - \
2*x[1]*Power(x[2],2)*r_X) - \
2*Power(a,2)*x[2]*Power(r,6)*(3*a*x[1]*r_Y + x[0]*(x[1] + 3*a*r_X)) + \
4*a*Power(r,8)*((Power(x[0],2) - Power(x[1],2))*r_Z - x[2]*(x[1]*r_Y \
+ x[0]*r_X)) + Power(a,6)*Power(x[2],2)*r*(-3*x[0]*x[1]*r_Z + \
x[2]*(2*a + 3*x[0]*r_Y - 3*x[1]*r_X)) + \
2*Power(a,2)*Power(r,7)*(3*x[0]*x[1]*r_Z + 2*x[2]*(a - x[0]*r_Y + \
x[1]*r_X)) + Power(r,9)*(-3*x[0]*x[1]*r_Z + x[2]*(2*a - 3*x[0]*r_Y + \
3*x[1]*r_X)) + Power(a,2)*Power(r,5)*(x[0]*x[1]*(Power(a,2) + \
Power(x[2],2))*r_Z + x[2]*(2*a*(Power(a,2) + Power(x[0],2) - \
Power(x[1],2) + Power(x[2],2)) + x[0]*(-Power(a,2) + \
Power(x[2],2))*r_Y + x[1]*(Power(a,2) - Power(x[2],2))*r_X)))) + \
2*x[4]*x[6]*r*(2*G*M*Power(r,4)*Power(a*x[1] + x[0]*r,2)*(a*x[0] - \
x[1]*r)*(4*Power(r,2)*(a*x[1] + x[0]*r)*(Power(a,2)*Power(x[2],2) + \
Power(r,4))*r_Z - 2*x[0]*r*(Power(a,2) + \
Power(r,2))*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_Z - 3*(a*x[1] + \
x[0]*r)*(Power(a,2) + Power(r,2))*(Power(a,2)*Power(x[2],2) + \
Power(r,4))*r_Z + 2*r*(a*x[1] + x[0]*r)*(Power(a,2) + \
Power(r,2))*(Power(a,2)*x[2] + 2*Power(r,3)*r_Z)) + \
2*G*M*Power(x[2],3)*r*(-(a*x[0]) + x[1]*r)*Power(Power(a,2) + \
Power(r,2),4)*(Power(a,2)*Power(x[2],2) - 3*Power(r,4))*r_X - \
2*G*M*Power(r,3)*(a*x[0] - x[1]*r)*Power(Power(a,2) + \
Power(r,2),2)*(2*Power(r,3)*(a*x[1] + \
x[0]*r)*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_Z - \
x[0]*Power(r,2)*(Power(a,2) + Power(r,2))*(Power(a,2)*Power(x[2],2) + \
Power(r,4))*r_Z - 3*r*(a*x[1] + x[0]*r)*(Power(a,2) + \
Power(r,2))*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_Z + \
2*Power(r,2)*(a*x[1] + x[0]*r)*(Power(a,2) + \
Power(r,2))*(Power(a,2)*x[2] + 2*Power(r,3)*r_Z) + \
4*x[2]*Power(r,4)*Power(Power(a,2) + Power(r,2),2)*r_X - \
2*x[2]*Power(Power(a,2) + Power(r,2),2)*(Power(a,2)*Power(x[2],2) + \
Power(r,4))*r_X) + (2*G*M*Power(r,3)*Power(a*x[0] - x[1]*r,2) - \
Power(Power(a,2) + Power(r,2),2)*(Power(a,2)*Power(x[2],2) + \
Power(r,4)))*(-2*Power(a,7)*Power(x[2],3)*(x[1]*r_Y + x[0]*r_X) + \
2*Power(a,5)*x[2]*Power(r,4)*(x[1]*r_Y + x[0]*r_X) + \
2*Power(a,2)*x[2]*Power(r,6)*(-(x[0]*x[1]) + 3*a*x[1]*r_Y + \
3*a*x[0]*r_X) + 2*Power(a,5)*x[2]*Power(r,2)*(a*x[0]*x[1] + \
2*(-Power(x[0],2) + Power(x[1],2))*x[2]*r_Z - x[1]*Power(x[2],2)*r_Y \
- x[0]*Power(x[2],2)*r_X) + \
2*Power(a,4)*x[2]*Power(r,3)*(a*Power(x[0],2) - a*Power(x[1],2) - \
2*a*Power(x[2],2) + 3*x[0]*x[1]*x[2]*r_Z - 2*x[0]*Power(x[2],2)*r_Y + \
2*x[1]*Power(x[2],2)*r_X) + 4*a*Power(r,8)*((Power(x[0],2) - \
Power(x[1],2))*r_Z + x[2]*(x[1]*r_Y + x[0]*r_X)) - \
Power(a,6)*Power(x[2],2)*r*(3*x[0]*x[1]*r_Z + x[2]*(2*a + 3*x[0]*r_Y \
- 3*x[1]*r_X)) + 2*Power(a,2)*Power(r,7)*(3*x[0]*x[1]*r_Z - 2*x[2]*(a \
- x[0]*r_Y + x[1]*r_X)) - Power(r,9)*(3*x[0]*x[1]*r_Z + x[2]*(2*a - \
3*x[0]*r_Y + 3*x[1]*r_X)) + \
Power(a,2)*Power(r,5)*(x[0]*x[1]*(Power(a,2) + Power(x[2],2))*r_Z - \
x[2]*(2*a*(Power(a,2) - Power(x[0],2) + Power(x[1],2) + \
Power(x[2],2)) + x[0]*(-Power(a,2) + Power(x[2],2))*r_Y + \
x[1]*(Power(a,2) - Power(x[2],2))*r_X))))))/(Power(Power(a,2) + \
Power(r,2),5)*Power(Power(a,2)*Power(x[2],2) + Power(r,4),3));
        dxdt[6] = (G*M*(Power(x[7],2)*Power(r,2)*Power(Power(a,2) + \
Power(r,2),3)*(-((Power(a,2) + Power(r,2))*(Power(a,2)*Power(x[2],2) \
- 2*G*M*Power(x[2],2)*r + Power(r,4))*(2*Power(a,2)*x[2]*r - \
3*Power(a,2)*Power(x[2],2)*r_Z + Power(r,4)*r_Z)) + \
2*G*M*x[2]*(a*x[0] - x[1]*r)*(3*Power(a,2)*Power(x[2],2)*Power(r,2) - \
Power(r,6))*r_Y - 2*G*M*x[2]*(a*x[1] + \
x[0]*r)*(3*Power(a,2)*Power(x[2],2)*Power(r,2) - Power(r,6))*r_X) - \
4*G*M*x[2]*x[7]*x[6]*Power(r,3)*(Power(a,2) + \
Power(r,2))*(-(r*Power(Power(a,2) + \
Power(r,2),3)*(2*Power(a,2)*x[2]*r - 3*Power(a,2)*Power(x[2],2)*r_Z + \
Power(r,4)*r_Z)) + (a*x[0] - x[1]*r)*(2*Power(r,3)*(a*x[0] - \
x[1]*r)*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_Z + \
x[1]*Power(r,2)*(Power(a,2) + Power(r,2))*(Power(a,2)*Power(x[2],2) + \
Power(r,4))*r_Z + 3*r*(-(a*x[0]) + x[1]*r)*(Power(a,2) + \
Power(r,2))*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_Z + \
Power(r,2)*(a*x[0] - x[1]*r)*(Power(a,2) + \
Power(r,2))*(2*Power(a,2)*x[2] + 4*Power(r,3)*r_Z) + \
4*x[2]*Power(r,4)*Power(Power(a,2) + Power(r,2),2)*r_Y - \
2*x[2]*Power(Power(a,2) + Power(r,2),2)*(Power(a,2)*Power(x[2],2) + \
Power(r,4))*r_Y) + (a*x[1] + x[0]*r)*(2*Power(r,3)*(a*x[1] + \
x[0]*r)*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_Z - \
x[0]*Power(r,2)*(Power(a,2) + Power(r,2))*(Power(a,2)*Power(x[2],2) + \
Power(r,4))*r_Z - 3*r*(a*x[1] + x[0]*r)*(Power(a,2) + \
Power(r,2))*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_Z + \
2*Power(r,2)*(a*x[1] + x[0]*r)*(Power(a,2) + \
Power(r,2))*(Power(a,2)*x[2] + 2*Power(r,3)*r_Z) - \
4*x[2]*Power(r,4)*Power(Power(a,2) + Power(r,2),2)*r_X + \
2*x[2]*Power(Power(a,2) + Power(r,2),2)*(Power(a,2)*Power(x[2],2) + \
Power(r,4))*r_X)) - x[2]*Power(x[6],2)*(Power(a,2) + \
Power(r,2))*(4*G*M*Power(r,3)*Power(Power(a,2) + \
Power(r,2),3)*(-(Power(a,2)*Power(x[2],2)) + Power(r,4))*(r - \
2*x[2]*r_Z) + Power(Power(a,2) + \
Power(r,2),3)*(Power(a,2)*Power(x[2],2) - 2*G*M*Power(x[2],2)*r + \
Power(r,4))*(2*Power(r,5) + Power(a,2)*Power(x[2],3)*r_Z - \
3*x[2]*Power(r,4)*r_Z) - 2*G*M*Power(r,2)*(a*x[0] - \
x[1]*r)*(2*Power(a,4)*x[1]*Power(x[2],2)*Power(r,3) - \
2*x[1]*Power(r,9) + 4*Power(a,5)*x[0]*Power(x[2],3)*r*r_Z + \
2*Power(a,2)*x[2]*Power(r,5)*(x[1]*x[2] - 2*a*x[0]*r_Z) - \
2*a*Power(r,7)*(a*x[1] + 4*x[0]*x[2]*r_Z) + \
Power(a,6)*Power(x[2],4)*r_Y + 2*Power(a,2)*Power(r,6)*(a*x[0] + \
x[1]*x[2]*r_Z - 3*Power(x[2],2)*r_Y) + Power(r,8)*(2*a*x[0] + \
6*x[1]*x[2]*r_Z - 3*Power(x[2],2)*r_Y) - \
2*Power(a,4)*Power(x[2],2)*Power(r,2)*(a*x[0] + 3*x[1]*x[2]*r_Z - \
Power(x[2],2)*r_Y) + Power(a,2)*Power(x[2],2)*Power(r,4)*(-2*a*x[0] - \
2*x[1]*x[2]*r_Z + (-3*Power(a,2) + Power(x[2],2))*r_Y)) - \
2*G*M*Power(r,2)*(a*x[1] + x[0]*r)*(2*Power(r,2)*(a*x[1] + \
x[0]*r)*(Power(a,2) + Power(r,2))*(Power(a,2)*Power(x[2],2) + \
Power(r,4)) - 4*x[2]*Power(r,3)*(a*x[1] + \
x[0]*r)*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_Z + \
2*x[0]*x[2]*Power(r,2)*(Power(a,2) + \
Power(r,2))*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_Z + \
4*x[2]*r*(a*x[1] + x[0]*r)*(Power(a,2) + \
Power(r,2))*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_Z - \
4*x[2]*Power(r,2)*(a*x[1] + x[0]*r)*(Power(a,2) + \
Power(r,2))*(Power(a,2)*x[2] + 2*Power(r,3)*r_Z) + \
4*Power(x[2],2)*Power(r,4)*Power(Power(a,2) + Power(r,2),2)*r_X - \
Power(x[2],2)*Power(Power(a,2) + \
Power(r,2),2)*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_X)) + \
2*x[7]*x[5]*r*(Power(a,2) + \
Power(r,2))*(-2*G*M*x[2]*r*Power(Power(a,2) + \
Power(r,2),3)*(3*Power(a,2)*Power(x[2],2)*Power(r,2) - \
Power(r,6))*r_Y + (Power(a,2) + Power(r,2))*(Power(a,2)*Power(x[2],2) \
- 2*G*M*Power(x[2],2)*r + Power(r,4))*(2*Power(r,3)*(a*x[0] - \
x[1]*r)*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_Z + \
x[1]*Power(r,2)*(Power(a,2) + Power(r,2))*(Power(a,2)*Power(x[2],2) + \
Power(r,4))*r_Z + 3*r*(-(a*x[0]) + x[1]*r)*(Power(a,2) + \
Power(r,2))*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_Z + \
Power(r,2)*(a*x[0] - x[1]*r)*(Power(a,2) + \
Power(r,2))*(2*Power(a,2)*x[2] + 4*Power(r,3)*r_Z) + \
4*x[2]*Power(r,4)*Power(Power(a,2) + Power(r,2),2)*r_Y - \
2*x[2]*Power(Power(a,2) + Power(r,2),2)*(Power(a,2)*Power(x[2],2) + \
Power(r,4))*r_Y) + 2*G*M*x[2]*Power(r,3)*(a*x[1] + \
x[0]*r)*(2*Power(a,3)*Power(r,5) + \
3*Power(a,5)*Power(x[2],2)*(x[1]*r_Y + x[0]*r_X) + \
Power(a,3)*Power(x[2],2)*Power(r,2)*(x[1]*r_Y + x[0]*r_X) - \
Power(a,3)*Power(r,4)*(x[1]*r_Y + x[0]*r_X) - \
3*a*Power(r,6)*(x[1]*r_Y + x[0]*r_X) + \
2*Power(a,4)*Power(x[2],2)*r*(a + 2*x[0]*r_Y - 2*x[1]*r_X) + \
2*Power(a,2)*Power(x[2],2)*Power(r,3)*(a + x[0]*r_Y - x[1]*r_X) + \
2*Power(r,7)*(a - x[0]*r_Y + x[1]*r_X))) + 2*x[7]*x[4]*r*(Power(a,2) \
+ Power(r,2))*(-2*G*M*x[2]*r*Power(Power(a,2) + \
Power(r,2),3)*(3*Power(a,2)*Power(x[2],2)*Power(r,2) - \
Power(r,6))*r_X - (Power(a,2) + Power(r,2))*(Power(a,2)*Power(x[2],2) \
- 2*G*M*Power(x[2],2)*r + Power(r,4))*(2*Power(r,3)*(a*x[1] + \
x[0]*r)*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_Z - \
x[0]*Power(r,2)*(Power(a,2) + Power(r,2))*(Power(a,2)*Power(x[2],2) + \
Power(r,4))*r_Z - 3*r*(a*x[1] + x[0]*r)*(Power(a,2) + \
Power(r,2))*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_Z + \
2*Power(r,2)*(a*x[1] + x[0]*r)*(Power(a,2) + \
Power(r,2))*(Power(a,2)*x[2] + 2*Power(r,3)*r_Z) - \
4*x[2]*Power(r,4)*Power(Power(a,2) + Power(r,2),2)*r_X + \
2*x[2]*Power(Power(a,2) + Power(r,2),2)*(Power(a,2)*Power(x[2],2) + \
Power(r,4))*r_X) + 2*G*M*x[2]*Power(r,3)*(a*x[0] - \
x[1]*r)*(2*Power(a,3)*Power(r,5) + \
3*Power(a,5)*Power(x[2],2)*(x[1]*r_Y + x[0]*r_X) + \
Power(a,3)*Power(x[2],2)*Power(r,2)*(x[1]*r_Y + x[0]*r_X) - \
Power(a,3)*Power(r,4)*(x[1]*r_Y + x[0]*r_X) - \
3*a*Power(r,6)*(x[1]*r_Y + x[0]*r_X) + \
2*Power(a,4)*Power(x[2],2)*r*(a + 2*x[0]*r_Y - 2*x[1]*r_X) + \
2*Power(a,2)*Power(x[2],2)*Power(r,3)*(a + x[0]*r_Y - x[1]*r_X) + \
2*Power(r,7)*(a - x[0]*r_Y + x[1]*r_X))) + \
Power(x[4],2)*r*(-2*G*M*x[2]*Power(r,3)*Power(a*x[1] + \
x[0]*r,2)*(4*Power(r,4)*(a*x[1] + x[0]*r)*(Power(a,2) + \
Power(r,2))*r_X + 4*Power(r,2)*(a*x[1] + \
x[0]*r)*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_X - 3*(a*x[1] + \
x[0]*r)*(Power(a,2) + Power(r,2))*(Power(a,2)*Power(x[2],2) + \
Power(r,4))*r_X - 2*r*(Power(a,2) + \
Power(r,2))*(Power(a,2)*Power(x[2],2) + Power(r,4))*(r + x[0]*r_X)) + \
4*G*M*x[2]*Power(r,3)*Power(Power(a,2) + \
Power(r,2),2)*(4*Power(r,4)*(a*x[1] + x[0]*r)*(Power(a,2) + \
Power(r,2))*r_X + 2*Power(r,2)*(a*x[1] + \
x[0]*r)*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_X - 3*(a*x[1] + \
x[0]*r)*(Power(a,2) + Power(r,2))*(Power(a,2)*Power(x[2],2) + \
Power(r,4))*r_X - r*(Power(a,2) + \
Power(r,2))*(Power(a,2)*Power(x[2],2) + Power(r,4))*(r + x[0]*r_X)) - \
(Power(a,2) + Power(r,2))*(Power(a,2)*Power(x[2],2) - \
2*G*M*Power(x[2],2)*r + Power(r,4))*(4*Power(r,3)*Power(a*x[1] + \
x[0]*r,2)*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_Z - \
2*x[0]*Power(r,2)*(a*x[1] + x[0]*r)*(Power(a,2) + \
Power(r,2))*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_Z - \
3*r*Power(a*x[1] + x[0]*r,2)*(Power(a,2) + \
Power(r,2))*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_Z + \
2*Power(r,2)*Power(a*x[1] + x[0]*r,2)*(Power(a,2) + \
Power(r,2))*(Power(a,2)*x[2] + 2*Power(r,3)*r_Z) - \
8*x[2]*Power(r,4)*(a*x[1] + x[0]*r)*Power(Power(a,2) + \
Power(r,2),2)*r_X - 4*x[2]*Power(r,2)*(a*x[1] + x[0]*r)*(Power(a,2) + \
Power(r,2))*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_X + \
4*x[2]*(a*x[1] + x[0]*r)*Power(Power(a,2) + \
Power(r,2),2)*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_X + \
2*x[2]*r*Power(Power(a,2) + Power(r,2),2)*(Power(a,2)*Power(x[2],2) + \
Power(r,4))*(r + x[0]*r_X)) - 2*G*M*x[2]*Power(r,3)*(a*x[0] - \
x[1]*r)*(-2*Power(a,4)*x[1]*Power(x[2],2)*Power(r,3) + \
2*Power(a,2)*x[1]*(-2*Power(a,2) + Power(x[2],2))*Power(r,5) + \
2*x[1]*Power(r,9) - 3*Power(a,6)*x[1]*Power(x[2],2)*(x[1]*r_Y + \
2*x[0]*r_X) - 2*a*Power(r,7)*(a*x[1] - 4*x[0]*x[1]*r_Y - \
4*(Power(x[0],2) - Power(x[1],2))*r_X) - \
4*Power(a,5)*Power(x[2],2)*r*(a*x[1] + 2*x[0]*x[1]*r_Y + \
2*(Power(x[0],2) - Power(x[1],2))*r_X) + \
Power(a,2)*Power(r,4)*(-6*a*x[0]*Power(x[2],2) + \
(Power(a,2)*Power(x[1],2) - Power(x[0],2)*Power(x[2],2))*r_Y + \
2*x[0]*x[1]*(Power(a,2) + Power(x[2],2))*r_X) + \
Power(a,2)*Power(r,6)*(-((Power(x[0],2) - 5*Power(x[1],2))*r_Y) - \
6*x[0]*(a - 2*x[1]*r_X)) + \
Power(a,4)*Power(x[2],2)*Power(r,2)*((-5*Power(x[0],2) + \
Power(x[1],2))*r_Y - 6*x[0]*(a - 2*x[1]*r_X)) + \
3*x[0]*Power(r,8)*(x[0]*r_Y - 2*(a + x[1]*r_X)))) - \
Power(x[5],2)*r*(4*G*M*x[2]*Power(r,3)*Power(Power(a,2) + \
Power(r,2),2)*(Power(r,8) - 3*Power(a,5)*x[0]*Power(x[2],2)*r_Y + \
4*Power(a,4)*x[1]*Power(x[2],2)*r*r_Y + \
2*Power(a,2)*x[1]*Power(x[2],2)*Power(r,3)*r_Y - \
2*x[1]*Power(r,7)*r_Y + Power(a,3)*Power(x[2],2)*Power(r,2)*(a - \
x[0]*r_Y) + a*Power(r,6)*(a + 3*x[0]*r_Y) + \
Power(a,2)*Power(r,4)*(Power(x[2],2) + a*x[0]*r_Y)) - \
2*G*M*x[2]*Power(r,3)*Power(a*x[0] - x[1]*r,2)*(4*Power(r,4)*(a*x[0] \
- x[1]*r)*(Power(a,2) + Power(r,2))*r_Y + 4*Power(r,2)*(a*x[0] - \
x[1]*r)*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_Y - 3*(a*x[0] - \
x[1]*r)*(Power(a,2) + Power(r,2))*(Power(a,2)*Power(x[2],2) + \
Power(r,4))*r_Y + 2*r*(Power(a,2) + \
Power(r,2))*(Power(a,2)*Power(x[2],2) + Power(r,4))*(r + x[1]*r_Y)) + \
(Power(a,2) + Power(r,2))*(Power(a,2)*Power(x[2],2) - \
2*G*M*Power(x[2],2)*r + Power(r,4))*(4*Power(r,3)*Power(a*x[0] - \
x[1]*r,2)*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_Z + \
2*x[1]*Power(r,2)*(a*x[0] - x[1]*r)*(Power(a,2) + \
Power(r,2))*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_Z - \
3*r*Power(a*x[0] - x[1]*r,2)*(Power(a,2) + \
Power(r,2))*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_Z + \
2*Power(r,2)*Power(a*x[0] - x[1]*r,2)*(Power(a,2) + \
Power(r,2))*(Power(a,2)*x[2] + 2*Power(r,3)*r_Z) + \
8*x[2]*Power(r,4)*(a*x[0] - x[1]*r)*Power(Power(a,2) + \
Power(r,2),2)*r_Y + 4*x[2]*Power(r,2)*(a*x[0] - x[1]*r)*(Power(a,2) + \
Power(r,2))*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_Y - \
4*x[2]*(a*x[0] - x[1]*r)*Power(Power(a,2) + \
Power(r,2),2)*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_Y + \
2*x[2]*r*Power(Power(a,2) + Power(r,2),2)*(Power(a,2)*Power(x[2],2) + \
Power(r,4))*(r + x[1]*r_Y)) - 2*G*M*x[2]*Power(r,3)*(a*x[1] + \
x[0]*r)*(-2*Power(a,4)*x[0]*Power(x[2],2)*Power(r,3) + \
2*Power(a,2)*x[0]*(-2*Power(a,2) + Power(x[2],2))*Power(r,5) + \
2*x[0]*Power(r,9) - 3*Power(a,6)*x[0]*Power(x[2],2)*(2*x[1]*r_Y + \
x[0]*r_X) + 3*x[1]*Power(r,8)*(2*a - 2*x[0]*r_Y + x[1]*r_X) + \
Power(a,4)*Power(x[2],2)*Power(r,2)*(6*a*x[1] + 12*x[0]*x[1]*r_Y + \
(Power(x[0],2) - 5*Power(x[1],2))*r_X) + \
Power(a,2)*Power(r,6)*(6*a*x[1] + 12*x[0]*x[1]*r_Y + (5*Power(x[0],2) \
- Power(x[1],2))*r_X) + Power(a,2)*Power(r,4)*(6*a*x[1]*Power(x[2],2) \
+ 2*x[0]*x[1]*(Power(a,2) + Power(x[2],2))*r_Y + \
(Power(a,2)*Power(x[0],2) - Power(x[1],2)*Power(x[2],2))*r_X) - \
4*Power(a,5)*Power(x[2],2)*r*(2*(Power(x[0],2) - Power(x[1],2))*r_Y + \
x[0]*(a - 2*x[1]*r_X)) - 2*a*Power(r,7)*(-4*(Power(x[0],2) - \
Power(x[1],2))*r_Y + x[0]*(a + 4*x[1]*r_X)))) - \
2*x[4]*x[5]*r*(2*G*M*x[2]*Power(r,3)*Power(a*x[1] + \
x[0]*r,2)*(4*Power(r,4)*(a*x[1] + x[0]*r)*(Power(a,2) + \
Power(r,2))*r_Y + 4*Power(r,2)*(a*x[1] + \
x[0]*r)*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_Y - 3*(a*x[1] + \
x[0]*r)*(Power(a,2) + Power(r,2))*(Power(a,2)*Power(x[2],2) + \
Power(r,4))*r_Y - 2*r*(Power(a,2) + \
Power(r,2))*(Power(a,2)*Power(x[2],2) + Power(r,4))*(a + x[0]*r_Y)) + \
2*G*M*x[2]*Power(r,3)*Power(Power(a,2) + \
Power(r,2),2)*((3*Power(a,5)*x[1]*Power(x[2],2) + \
4*Power(a,4)*x[0]*Power(x[2],2)*r + \
Power(a,3)*x[1]*Power(x[2],2)*Power(r,2) + \
2*Power(a,2)*x[0]*Power(x[2],2)*Power(r,3) - \
Power(a,3)*x[1]*Power(r,4) - 3*a*x[1]*Power(r,6) - \
2*x[0]*Power(r,7))*r_Y + (-3*Power(a,5)*x[0]*Power(x[2],2) + \
4*Power(a,4)*x[1]*Power(x[2],2)*r - \
Power(a,3)*x[0]*Power(x[2],2)*Power(r,2) + \
2*Power(a,2)*x[1]*Power(x[2],2)*Power(r,3) + \
Power(a,3)*x[0]*Power(r,4) + 3*a*x[0]*Power(r,6) - \
2*x[1]*Power(r,7))*r_X) - 2*G*M*x[2]*Power(r,3)*Power(a*x[0] - \
x[1]*r,2)*(4*Power(r,4)*(a*x[0] - x[1]*r)*(Power(a,2) + \
Power(r,2))*r_X + 4*Power(r,2)*(a*x[0] - \
x[1]*r)*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_X - 3*(a*x[0] - \
x[1]*r)*(Power(a,2) + Power(r,2))*(Power(a,2)*Power(x[2],2) + \
Power(r,4))*r_X - 2*r*(Power(a,2) + \
Power(r,2))*(Power(a,2)*Power(x[2],2) + Power(r,4))*(a - x[1]*r_X)) - \
(Power(a,2) + Power(r,2))*(Power(a,2)*Power(x[2],2) - \
2*G*M*Power(x[2],2)*r + \
Power(r,4))*(2*Power(a,5)*x[2]*Power(r,4)*(x[1]*r_Y - x[0]*r_X) + \
2*Power(a,7)*Power(x[2],3)*(-(x[1]*r_Y) + x[0]*r_X) - \
3*Power(a,6)*Power(x[2],2)*r*(x[0]*x[1]*r_Z + x[0]*x[2]*r_Y + \
x[1]*x[2]*r_X) + 2*Power(a,5)*x[2]*Power(r,2)*(a*x[0]*x[1] + \
2*(-Power(x[0],2) + Power(x[1],2))*x[2]*r_Z - x[1]*Power(x[2],2)*r_Y \
+ x[0]*Power(x[2],2)*r_X) + \
2*Power(a,4)*x[2]*Power(r,3)*(a*Power(x[0],2) - a*Power(x[1],2) + \
3*x[0]*x[1]*x[2]*r_Z - 2*x[0]*Power(x[2],2)*r_Y - \
2*x[1]*Power(x[2],2)*r_X) - \
2*Power(a,2)*x[2]*Power(r,6)*(-3*a*x[1]*r_Y + x[0]*(x[1] + 3*a*r_X)) \
+ 4*a*Power(r,8)*((Power(x[0],2) - Power(x[1],2))*r_Z + \
x[2]*(x[1]*r_Y - x[0]*r_X)) + \
2*Power(a,2)*Power(r,7)*(3*x[0]*x[1]*r_Z + 2*x[2]*(x[0]*r_Y + \
x[1]*r_X)) + Power(r,9)*(-3*x[0]*x[1]*r_Z + 3*x[2]*(x[0]*r_Y + \
x[1]*r_X)) + Power(a,2)*Power(r,5)*(x[0]*x[1]*(Power(a,2) + \
Power(x[2],2))*r_Z + x[2]*(2*a*(Power(x[0],2) - Power(x[1],2)) + \
x[0]*(Power(a,2) - Power(x[2],2))*r_Y + x[1]*(Power(a,2) - \
Power(x[2],2))*r_X)))) + \
2*x[2]*x[5]*x[6]*(2*G*M*Power(r,4)*Power(a*x[0] - \
x[1]*r,2)*(4*Power(r,2)*(a*x[0] - x[1]*r)*(Power(a,2)*Power(x[2],2) + \
Power(r,4))*r_Z + 2*x[1]*r*(Power(a,2) + \
Power(r,2))*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_Z - 3*(a*x[0] - \
x[1]*r)*(Power(a,2) + Power(r,2))*(Power(a,2)*Power(x[2],2) + \
Power(r,4))*r_Z + r*(a*x[0] - x[1]*r)*(Power(a,2) + \
Power(r,2))*(2*Power(a,2)*x[2] + 4*Power(r,3)*r_Z)) - \
x[2]*Power(Power(a,2) + Power(r,2),4)*(Power(a,2)*Power(x[2],2) - \
3*Power(r,4))*(Power(a,2)*Power(x[2],2) - 2*G*M*Power(x[2],2)*r + \
Power(r,4))*r_Y - 2*G*M*Power(r,3)*Power(Power(a,2) + \
Power(r,2),2)*(2*Power(r,3)*(a*x[0] - \
x[1]*r)*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_Z + \
x[1]*Power(r,2)*(Power(a,2) + Power(r,2))*(Power(a,2)*Power(x[2],2) + \
Power(r,4))*r_Z + 3*r*(-(a*x[0]) + x[1]*r)*(Power(a,2) + \
Power(r,2))*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_Z + \
Power(r,2)*(a*x[0] - x[1]*r)*(Power(a,2) + \
Power(r,2))*(2*Power(a,2)*x[2] + 4*Power(r,3)*r_Z) - \
4*x[2]*Power(r,4)*Power(Power(a,2) + Power(r,2),2)*r_Y + \
2*x[2]*Power(Power(a,2) + Power(r,2),2)*(Power(a,2)*Power(x[2],2) + \
Power(r,4))*r_Y) + 2*G*M*Power(r,3)*(a*x[1] + \
x[0]*r)*(2*Power(a,7)*Power(x[2],3)*(x[1]*r_Y + x[0]*r_X) - \
2*Power(a,5)*x[2]*Power(r,4)*(x[1]*r_Y + x[0]*r_X) + \
2*Power(a,5)*x[2]*Power(r,2)*(a*x[0]*x[1] + 2*(-Power(x[0],2) + \
Power(x[1],2))*x[2]*r_Z + x[1]*Power(x[2],2)*r_Y + \
x[0]*Power(x[2],2)*r_X) + \
2*Power(a,4)*x[2]*Power(r,3)*(a*Power(x[0],2) - a*Power(x[1],2) + \
2*a*Power(x[2],2) + 3*x[0]*x[1]*x[2]*r_Z + 2*x[0]*Power(x[2],2)*r_Y - \
2*x[1]*Power(x[2],2)*r_X) - \
2*Power(a,2)*x[2]*Power(r,6)*(3*a*x[1]*r_Y + x[0]*(x[1] + 3*a*r_X)) + \
4*a*Power(r,8)*((Power(x[0],2) - Power(x[1],2))*r_Z - x[2]*(x[1]*r_Y \
+ x[0]*r_X)) + Power(a,6)*Power(x[2],2)*r*(-3*x[0]*x[1]*r_Z + \
x[2]*(2*a + 3*x[0]*r_Y - 3*x[1]*r_X)) + \
2*Power(a,2)*Power(r,7)*(3*x[0]*x[1]*r_Z + 2*x[2]*(a - x[0]*r_Y + \
x[1]*r_X)) + Power(r,9)*(-3*x[0]*x[1]*r_Z + x[2]*(2*a - 3*x[0]*r_Y + \
3*x[1]*r_X)) + Power(a,2)*Power(r,5)*(x[0]*x[1]*(Power(a,2) + \
Power(x[2],2))*r_Z + x[2]*(2*a*(Power(a,2) + Power(x[0],2) - \
Power(x[1],2) + Power(x[2],2)) + x[0]*(-Power(a,2) + \
Power(x[2],2))*r_Y + x[1]*(Power(a,2) - Power(x[2],2))*r_X)))) - \
2*x[2]*x[4]*x[6]*(2*G*M*Power(r,4)*Power(a*x[1] + \
x[0]*r,2)*(4*Power(r,2)*(a*x[1] + x[0]*r)*(Power(a,2)*Power(x[2],2) + \
Power(r,4))*r_Z - 2*x[0]*r*(Power(a,2) + \
Power(r,2))*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_Z - 3*(a*x[1] + \
x[0]*r)*(Power(a,2) + Power(r,2))*(Power(a,2)*Power(x[2],2) + \
Power(r,4))*r_Z + 2*r*(a*x[1] + x[0]*r)*(Power(a,2) + \
Power(r,2))*(Power(a,2)*x[2] + 2*Power(r,3)*r_Z)) + \
x[2]*Power(Power(a,2) + Power(r,2),4)*(Power(a,2)*Power(x[2],2) - \
3*Power(r,4))*(Power(a,2)*Power(x[2],2) - 2*G*M*Power(x[2],2)*r + \
Power(r,4))*r_X - 2*G*M*Power(r,3)*Power(Power(a,2) + \
Power(r,2),2)*(2*Power(r,3)*(a*x[1] + \
x[0]*r)*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_Z - \
x[0]*Power(r,2)*(Power(a,2) + Power(r,2))*(Power(a,2)*Power(x[2],2) + \
Power(r,4))*r_Z - 3*r*(a*x[1] + x[0]*r)*(Power(a,2) + \
Power(r,2))*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_Z + \
2*Power(r,2)*(a*x[1] + x[0]*r)*(Power(a,2) + \
Power(r,2))*(Power(a,2)*x[2] + 2*Power(r,3)*r_Z) + \
4*x[2]*Power(r,4)*Power(Power(a,2) + Power(r,2),2)*r_X - \
2*x[2]*Power(Power(a,2) + Power(r,2),2)*(Power(a,2)*Power(x[2],2) + \
Power(r,4))*r_X) + 2*G*M*Power(r,3)*(a*x[0] - \
x[1]*r)*(-2*Power(a,7)*Power(x[2],3)*(x[1]*r_Y + x[0]*r_X) + \
2*Power(a,5)*x[2]*Power(r,4)*(x[1]*r_Y + x[0]*r_X) + \
2*Power(a,2)*x[2]*Power(r,6)*(-(x[0]*x[1]) + 3*a*x[1]*r_Y + \
3*a*x[0]*r_X) + 2*Power(a,5)*x[2]*Power(r,2)*(a*x[0]*x[1] + \
2*(-Power(x[0],2) + Power(x[1],2))*x[2]*r_Z - x[1]*Power(x[2],2)*r_Y \
- x[0]*Power(x[2],2)*r_X) + \
2*Power(a,4)*x[2]*Power(r,3)*(a*Power(x[0],2) - a*Power(x[1],2) - \
2*a*Power(x[2],2) + 3*x[0]*x[1]*x[2]*r_Z - 2*x[0]*Power(x[2],2)*r_Y + \
2*x[1]*Power(x[2],2)*r_X) + 4*a*Power(r,8)*((Power(x[0],2) - \
Power(x[1],2))*r_Z + x[2]*(x[1]*r_Y + x[0]*r_X)) - \
Power(a,6)*Power(x[2],2)*r*(3*x[0]*x[1]*r_Z + x[2]*(2*a + 3*x[0]*r_Y \
- 3*x[1]*r_X)) + 2*Power(a,2)*Power(r,7)*(3*x[0]*x[1]*r_Z - 2*x[2]*(a \
- x[0]*r_Y + x[1]*r_X)) - Power(r,9)*(3*x[0]*x[1]*r_Z + x[2]*(2*a - \
3*x[0]*r_Y + 3*x[1]*r_X)) + \
Power(a,2)*Power(r,5)*(x[0]*x[1]*(Power(a,2) + Power(x[2],2))*r_Z - \
x[2]*(2*a*(Power(a,2) - Power(x[0],2) + Power(x[1],2) + \
Power(x[2],2)) + x[0]*(-Power(a,2) + Power(x[2],2))*r_Y + \
x[1]*(Power(a,2) - Power(x[2],2))*r_X))))))/(Power(Power(a,2) + \
Power(r,2),4)*Power(Power(a,2)*Power(x[2],2) + Power(r,4),3));
        dxdt[7] = (2*G*M*(G*M*Power(x[7],2)*Power(r,3)*Power(Power(a,2) + \
Power(r,2),3)*(x[2]*r*(Power(a,2) + Power(r,2))*(-2*Power(a,2)*x[2]*r \
+ 3*Power(a,2)*Power(x[2],2)*r_Z - Power(r,4)*r_Z) - (-(a*x[0]) + \
x[1]*r)*(-3*Power(a,2)*Power(x[2],2)*Power(r,2) + Power(r,6))*r_Y + \
(a*x[1] + x[0]*r)*(3*Power(a,2)*Power(x[2],2)*Power(r,2) - \
Power(r,6))*r_X) + x[7]*x[6]*(Power(a,2) + \
Power(r,2))*(-(Power(r,2)*Power(Power(a,2) + \
Power(r,2),3)*(Power(a,2)*Power(x[2],2) + 2*G*M*Power(r,3) + \
Power(r,4))*(2*Power(a,2)*x[2]*r - 3*Power(a,2)*Power(x[2],2)*r_Z + \
Power(r,4)*r_Z)) + 2*G*M*Power(r,4)*(a*x[0] - \
x[1]*r)*(2*Power(r,3)*(a*x[0] - x[1]*r)*(Power(a,2)*Power(x[2],2) + \
Power(r,4))*r_Z + x[1]*Power(r,2)*(Power(a,2) + \
Power(r,2))*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_Z + \
3*r*(-(a*x[0]) + x[1]*r)*(Power(a,2) + \
Power(r,2))*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_Z + \
Power(r,2)*(a*x[0] - x[1]*r)*(Power(a,2) + \
Power(r,2))*(2*Power(a,2)*x[2] + 4*Power(r,3)*r_Z) + \
4*x[2]*Power(r,4)*Power(Power(a,2) + Power(r,2),2)*r_Y - \
2*x[2]*Power(Power(a,2) + Power(r,2),2)*(Power(a,2)*Power(x[2],2) + \
Power(r,4))*r_Y) + 2*G*M*Power(r,4)*(a*x[1] + \
x[0]*r)*(2*Power(r,3)*(a*x[1] + x[0]*r)*(Power(a,2)*Power(x[2],2) + \
Power(r,4))*r_Z - x[0]*Power(r,2)*(Power(a,2) + \
Power(r,2))*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_Z - 3*r*(a*x[1] \
+ x[0]*r)*(Power(a,2) + Power(r,2))*(Power(a,2)*Power(x[2],2) + \
Power(r,4))*r_Z + 2*Power(r,2)*(a*x[1] + x[0]*r)*(Power(a,2) + \
Power(r,2))*(Power(a,2)*x[2] + 2*Power(r,3)*r_Z) - \
4*x[2]*Power(r,4)*Power(Power(a,2) + Power(r,2),2)*r_X + \
2*x[2]*Power(Power(a,2) + Power(r,2),2)*(Power(a,2)*Power(x[2],2) + \
Power(r,4))*r_X)) - Power(x[6],2)*r*(Power(a,2) + \
Power(r,2))*(Power(Power(a,2) + \
Power(r,2),3)*(Power(a,2)*Power(x[2],2) - \
Power(r,4))*(Power(a,2)*Power(x[2],2) + 2*G*M*Power(r,3) + \
Power(r,4))*(r - 2*x[2]*r_Z) + G*M*Power(x[2],2)*r*Power(Power(a,2) + \
Power(r,2),3)*(2*Power(r,5) + Power(a,2)*Power(x[2],3)*r_Z - \
3*x[2]*Power(r,4)*r_Z) + G*M*Power(r,2)*(a*x[0] - \
x[1]*r)*(2*Power(a,4)*x[1]*Power(x[2],2)*Power(r,3) - \
2*x[1]*Power(r,9) + 4*Power(a,5)*x[0]*Power(x[2],3)*r*r_Z + \
2*Power(a,2)*x[2]*Power(r,5)*(x[1]*x[2] - 2*a*x[0]*r_Z) - \
2*a*Power(r,7)*(a*x[1] + 4*x[0]*x[2]*r_Z) + \
Power(a,6)*Power(x[2],4)*r_Y + 2*Power(a,2)*Power(r,6)*(a*x[0] + \
x[1]*x[2]*r_Z - 3*Power(x[2],2)*r_Y) + Power(r,8)*(2*a*x[0] + \
6*x[1]*x[2]*r_Z - 3*Power(x[2],2)*r_Y) - \
2*Power(a,4)*Power(x[2],2)*Power(r,2)*(a*x[0] + 3*x[1]*x[2]*r_Z - \
Power(x[2],2)*r_Y) + Power(a,2)*Power(x[2],2)*Power(r,4)*(-2*a*x[0] - \
2*x[1]*x[2]*r_Z + (-3*Power(a,2) + Power(x[2],2))*r_Y)) + \
G*M*Power(r,2)*(a*x[1] + x[0]*r)*(2*Power(r,2)*(a*x[1] + \
x[0]*r)*(Power(a,2) + Power(r,2))*(Power(a,2)*Power(x[2],2) + \
Power(r,4)) - 4*x[2]*Power(r,3)*(a*x[1] + \
x[0]*r)*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_Z + \
2*x[0]*x[2]*Power(r,2)*(Power(a,2) + \
Power(r,2))*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_Z + \
4*x[2]*r*(a*x[1] + x[0]*r)*(Power(a,2) + \
Power(r,2))*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_Z - \
4*x[2]*Power(r,2)*(a*x[1] + x[0]*r)*(Power(a,2) + \
Power(r,2))*(Power(a,2)*x[2] + 2*Power(r,3)*r_Z) + \
4*Power(x[2],2)*Power(r,4)*Power(Power(a,2) + Power(r,2),2)*r_X - \
Power(x[2],2)*Power(Power(a,2) + \
Power(r,2),2)*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_X)) - \
x[7]*x[5]*(Power(a,2) + Power(r,2))*(-(Power(Power(a,2) + \
Power(r,2),3)*(Power(a,2)*Power(x[2],2) + 2*G*M*Power(r,3) + \
Power(r,4))*(3*Power(a,2)*Power(x[2],2)*Power(r,2) - Power(r,6))*r_Y) \
- 2*G*M*x[2]*Power(r,3)*(Power(a,2) + \
Power(r,2))*(2*Power(r,3)*(a*x[0] - x[1]*r)*(Power(a,2)*Power(x[2],2) \
+ Power(r,4))*r_Z + x[1]*Power(r,2)*(Power(a,2) + \
Power(r,2))*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_Z + \
3*r*(-(a*x[0]) + x[1]*r)*(Power(a,2) + \
Power(r,2))*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_Z + \
Power(r,2)*(a*x[0] - x[1]*r)*(Power(a,2) + \
Power(r,2))*(2*Power(a,2)*x[2] + 4*Power(r,3)*r_Z) + \
4*x[2]*Power(r,4)*Power(Power(a,2) + Power(r,2),2)*r_Y - \
2*x[2]*Power(Power(a,2) + Power(r,2),2)*(Power(a,2)*Power(x[2],2) + \
Power(r,4))*r_Y) + 2*G*M*Power(r,5)*(a*x[1] + \
x[0]*r)*(2*Power(a,3)*Power(r,5) + \
3*Power(a,5)*Power(x[2],2)*(x[1]*r_Y + x[0]*r_X) + \
Power(a,3)*Power(x[2],2)*Power(r,2)*(x[1]*r_Y + x[0]*r_X) - \
Power(a,3)*Power(r,4)*(x[1]*r_Y + x[0]*r_X) - \
3*a*Power(r,6)*(x[1]*r_Y + x[0]*r_X) + \
2*Power(a,4)*Power(x[2],2)*r*(a + 2*x[0]*r_Y - 2*x[1]*r_X) + \
2*Power(a,2)*Power(x[2],2)*Power(r,3)*(a + x[0]*r_Y - x[1]*r_X) + \
2*Power(r,7)*(a - x[0]*r_Y + x[1]*r_X))) - x[7]*x[4]*(Power(a,2) + \
Power(r,2))*(-(Power(Power(a,2) + \
Power(r,2),3)*(Power(a,2)*Power(x[2],2) + 2*G*M*Power(r,3) + \
Power(r,4))*(3*Power(a,2)*Power(x[2],2)*Power(r,2) - Power(r,6))*r_X) \
+ 2*G*M*x[2]*Power(r,3)*(Power(a,2) + \
Power(r,2))*(2*Power(r,3)*(a*x[1] + x[0]*r)*(Power(a,2)*Power(x[2],2) \
+ Power(r,4))*r_Z - x[0]*Power(r,2)*(Power(a,2) + \
Power(r,2))*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_Z - 3*r*(a*x[1] \
+ x[0]*r)*(Power(a,2) + Power(r,2))*(Power(a,2)*Power(x[2],2) + \
Power(r,4))*r_Z + 2*Power(r,2)*(a*x[1] + x[0]*r)*(Power(a,2) + \
Power(r,2))*(Power(a,2)*x[2] + 2*Power(r,3)*r_Z) - \
4*x[2]*Power(r,4)*Power(Power(a,2) + Power(r,2),2)*r_X + \
2*x[2]*Power(Power(a,2) + Power(r,2),2)*(Power(a,2)*Power(x[2],2) + \
Power(r,4))*r_X) + 2*G*M*Power(r,5)*(a*x[0] - \
x[1]*r)*(2*Power(a,3)*Power(r,5) + \
3*Power(a,5)*Power(x[2],2)*(x[1]*r_Y + x[0]*r_X) + \
Power(a,3)*Power(x[2],2)*Power(r,2)*(x[1]*r_Y + x[0]*r_X) - \
Power(a,3)*Power(r,4)*(x[1]*r_Y + x[0]*r_X) - \
3*a*Power(r,6)*(x[1]*r_Y + x[0]*r_X) + \
2*Power(a,4)*Power(x[2],2)*r*(a + 2*x[0]*r_Y - 2*x[1]*r_X) + \
2*Power(a,2)*Power(x[2],2)*Power(r,3)*(a + x[0]*r_Y - x[1]*r_X) + \
2*Power(r,7)*(a - x[0]*r_Y + x[1]*r_X))) - \
Power(x[4],2)*Power(r,2)*(-(G*M*Power(r,3)*Power(a*x[1] + \
x[0]*r,2)*(4*Power(r,4)*(a*x[1] + x[0]*r)*(Power(a,2) + \
Power(r,2))*r_X + 4*Power(r,2)*(a*x[1] + \
x[0]*r)*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_X - 3*(a*x[1] + \
x[0]*r)*(Power(a,2) + Power(r,2))*(Power(a,2)*Power(x[2],2) + \
Power(r,4))*r_X - 2*r*(Power(a,2) + \
Power(r,2))*(Power(a,2)*Power(x[2],2) + Power(r,4))*(r + x[0]*r_X))) \
+ Power(Power(a,2) + Power(r,2),2)*(Power(a,2)*Power(x[2],2) + \
2*G*M*Power(r,3) + Power(r,4))*(4*Power(r,4)*(a*x[1] + \
x[0]*r)*(Power(a,2) + Power(r,2))*r_X + 2*Power(r,2)*(a*x[1] + \
x[0]*r)*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_X - 3*(a*x[1] + \
x[0]*r)*(Power(a,2) + Power(r,2))*(Power(a,2)*Power(x[2],2) + \
Power(r,4))*r_X - r*(Power(a,2) + \
Power(r,2))*(Power(a,2)*Power(x[2],2) + Power(r,4))*(r + x[0]*r_X)) + \
G*M*x[2]*r*(Power(a,2) + Power(r,2))*(4*Power(r,3)*Power(a*x[1] + \
x[0]*r,2)*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_Z - \
2*x[0]*Power(r,2)*(a*x[1] + x[0]*r)*(Power(a,2) + \
Power(r,2))*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_Z - \
3*r*Power(a*x[1] + x[0]*r,2)*(Power(a,2) + \
Power(r,2))*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_Z + \
2*Power(r,2)*Power(a*x[1] + x[0]*r,2)*(Power(a,2) + \
Power(r,2))*(Power(a,2)*x[2] + 2*Power(r,3)*r_Z) - \
8*x[2]*Power(r,4)*(a*x[1] + x[0]*r)*Power(Power(a,2) + \
Power(r,2),2)*r_X - 4*x[2]*Power(r,2)*(a*x[1] + x[0]*r)*(Power(a,2) + \
Power(r,2))*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_X + \
4*x[2]*(a*x[1] + x[0]*r)*Power(Power(a,2) + \
Power(r,2),2)*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_X + \
2*x[2]*r*Power(Power(a,2) + Power(r,2),2)*(Power(a,2)*Power(x[2],2) + \
Power(r,4))*(r + x[0]*r_X)) - G*M*Power(r,3)*(a*x[0] - \
x[1]*r)*(-2*Power(a,4)*x[1]*Power(x[2],2)*Power(r,3) + \
2*Power(a,2)*x[1]*(-2*Power(a,2) + Power(x[2],2))*Power(r,5) + \
2*x[1]*Power(r,9) - 3*Power(a,6)*x[1]*Power(x[2],2)*(x[1]*r_Y + \
2*x[0]*r_X) - 2*a*Power(r,7)*(a*x[1] - 4*x[0]*x[1]*r_Y - \
4*(Power(x[0],2) - Power(x[1],2))*r_X) - \
4*Power(a,5)*Power(x[2],2)*r*(a*x[1] + 2*x[0]*x[1]*r_Y + \
2*(Power(x[0],2) - Power(x[1],2))*r_X) + \
Power(a,2)*Power(r,4)*(-6*a*x[0]*Power(x[2],2) + \
(Power(a,2)*Power(x[1],2) - Power(x[0],2)*Power(x[2],2))*r_Y + \
2*x[0]*x[1]*(Power(a,2) + Power(x[2],2))*r_X) + \
Power(a,2)*Power(r,6)*(-((Power(x[0],2) - 5*Power(x[1],2))*r_Y) - \
6*x[0]*(a - 2*x[1]*r_X)) + \
Power(a,4)*Power(x[2],2)*Power(r,2)*((-5*Power(x[0],2) + \
Power(x[1],2))*r_Y - 6*x[0]*(a - 2*x[1]*r_X)) + \
3*x[0]*Power(r,8)*(x[0]*r_Y - 2*(a + x[1]*r_X)))) + \
Power(x[5],2)*Power(r,2)*(Power(Power(a,2) + \
Power(r,2),2)*(Power(a,2)*Power(x[2],2) + 2*G*M*Power(r,3) + \
Power(r,4))*(Power(r,8) - 3*Power(a,5)*x[0]*Power(x[2],2)*r_Y + \
4*Power(a,4)*x[1]*Power(x[2],2)*r*r_Y + \
2*Power(a,2)*x[1]*Power(x[2],2)*Power(r,3)*r_Y - \
2*x[1]*Power(r,7)*r_Y + Power(a,3)*Power(x[2],2)*Power(r,2)*(a - \
x[0]*r_Y) + a*Power(r,6)*(a + 3*x[0]*r_Y) + \
Power(a,2)*Power(r,4)*(Power(x[2],2) + a*x[0]*r_Y)) - \
G*M*Power(r,3)*Power(a*x[0] - x[1]*r,2)*(4*Power(r,4)*(a*x[0] - \
x[1]*r)*(Power(a,2) + Power(r,2))*r_Y + 4*Power(r,2)*(a*x[0] - \
x[1]*r)*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_Y - 3*(a*x[0] - \
x[1]*r)*(Power(a,2) + Power(r,2))*(Power(a,2)*Power(x[2],2) + \
Power(r,4))*r_Y + 2*r*(Power(a,2) + \
Power(r,2))*(Power(a,2)*Power(x[2],2) + Power(r,4))*(r + x[1]*r_Y)) - \
G*M*x[2]*r*(Power(a,2) + Power(r,2))*(4*Power(r,3)*Power(a*x[0] - \
x[1]*r,2)*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_Z + \
2*x[1]*Power(r,2)*(a*x[0] - x[1]*r)*(Power(a,2) + \
Power(r,2))*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_Z - \
3*r*Power(a*x[0] - x[1]*r,2)*(Power(a,2) + \
Power(r,2))*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_Z + \
2*Power(r,2)*Power(a*x[0] - x[1]*r,2)*(Power(a,2) + \
Power(r,2))*(Power(a,2)*x[2] + 2*Power(r,3)*r_Z) + \
8*x[2]*Power(r,4)*(a*x[0] - x[1]*r)*Power(Power(a,2) + \
Power(r,2),2)*r_Y + 4*x[2]*Power(r,2)*(a*x[0] - x[1]*r)*(Power(a,2) + \
Power(r,2))*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_Y - \
4*x[2]*(a*x[0] - x[1]*r)*Power(Power(a,2) + \
Power(r,2),2)*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_Y + \
2*x[2]*r*Power(Power(a,2) + Power(r,2),2)*(Power(a,2)*Power(x[2],2) + \
Power(r,4))*(r + x[1]*r_Y)) - G*M*Power(r,3)*(a*x[1] + \
x[0]*r)*(-2*Power(a,4)*x[0]*Power(x[2],2)*Power(r,3) + \
2*Power(a,2)*x[0]*(-2*Power(a,2) + Power(x[2],2))*Power(r,5) + \
2*x[0]*Power(r,9) - 3*Power(a,6)*x[0]*Power(x[2],2)*(2*x[1]*r_Y + \
x[0]*r_X) + 3*x[1]*Power(r,8)*(2*a - 2*x[0]*r_Y + x[1]*r_X) + \
Power(a,4)*Power(x[2],2)*Power(r,2)*(6*a*x[1] + 12*x[0]*x[1]*r_Y + \
(Power(x[0],2) - 5*Power(x[1],2))*r_X) + \
Power(a,2)*Power(r,6)*(6*a*x[1] + 12*x[0]*x[1]*r_Y + (5*Power(x[0],2) \
- Power(x[1],2))*r_X) + Power(a,2)*Power(r,4)*(6*a*x[1]*Power(x[2],2) \
+ 2*x[0]*x[1]*(Power(a,2) + Power(x[2],2))*r_Y + \
(Power(a,2)*Power(x[0],2) - Power(x[1],2)*Power(x[2],2))*r_X) - \
4*Power(a,5)*Power(x[2],2)*r*(2*(Power(x[0],2) - Power(x[1],2))*r_Y + \
x[0]*(a - 2*x[1]*r_X)) - 2*a*Power(r,7)*(-4*(Power(x[0],2) - \
Power(x[1],2))*r_Y + x[0]*(a + 4*x[1]*r_X)))) + \
x[4]*x[5]*Power(r,2)*(2*G*M*Power(r,3)*Power(a*x[1] + \
x[0]*r,2)*(4*Power(r,4)*(a*x[1] + x[0]*r)*(Power(a,2) + \
Power(r,2))*r_Y + 4*Power(r,2)*(a*x[1] + \
x[0]*r)*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_Y - 3*(a*x[1] + \
x[0]*r)*(Power(a,2) + Power(r,2))*(Power(a,2)*Power(x[2],2) + \
Power(r,4))*r_Y - 2*r*(Power(a,2) + \
Power(r,2))*(Power(a,2)*Power(x[2],2) + Power(r,4))*(a + x[0]*r_Y)) + \
Power(Power(a,2) + Power(r,2),2)*(Power(a,2)*Power(x[2],2) + \
2*G*M*Power(r,3) + Power(r,4))*((3*Power(a,5)*x[1]*Power(x[2],2) + \
4*Power(a,4)*x[0]*Power(x[2],2)*r + \
Power(a,3)*x[1]*Power(x[2],2)*Power(r,2) + \
2*Power(a,2)*x[0]*Power(x[2],2)*Power(r,3) - \
Power(a,3)*x[1]*Power(r,4) - 3*a*x[1]*Power(r,6) - \
2*x[0]*Power(r,7))*r_Y + (-3*Power(a,5)*x[0]*Power(x[2],2) + \
4*Power(a,4)*x[1]*Power(x[2],2)*r - \
Power(a,3)*x[0]*Power(x[2],2)*Power(r,2) + \
2*Power(a,2)*x[1]*Power(x[2],2)*Power(r,3) + \
Power(a,3)*x[0]*Power(r,4) + 3*a*x[0]*Power(r,6) - \
2*x[1]*Power(r,7))*r_X) - 2*G*M*Power(r,3)*Power(a*x[0] - \
x[1]*r,2)*(4*Power(r,4)*(a*x[0] - x[1]*r)*(Power(a,2) + \
Power(r,2))*r_X + 4*Power(r,2)*(a*x[0] - \
x[1]*r)*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_X - 3*(a*x[0] - \
x[1]*r)*(Power(a,2) + Power(r,2))*(Power(a,2)*Power(x[2],2) + \
Power(r,4))*r_X - 2*r*(Power(a,2) + \
Power(r,2))*(Power(a,2)*Power(x[2],2) + Power(r,4))*(a - x[1]*r_X)) + \
2*G*M*x[2]*r*(Power(a,2) + \
Power(r,2))*(2*Power(a,5)*x[2]*Power(r,4)*(x[1]*r_Y - x[0]*r_X) + \
2*Power(a,7)*Power(x[2],3)*(-(x[1]*r_Y) + x[0]*r_X) - \
3*Power(a,6)*Power(x[2],2)*r*(x[0]*x[1]*r_Z + x[0]*x[2]*r_Y + \
x[1]*x[2]*r_X) + 2*Power(a,5)*x[2]*Power(r,2)*(a*x[0]*x[1] + \
2*(-Power(x[0],2) + Power(x[1],2))*x[2]*r_Z - x[1]*Power(x[2],2)*r_Y \
+ x[0]*Power(x[2],2)*r_X) + \
2*Power(a,4)*x[2]*Power(r,3)*(a*Power(x[0],2) - a*Power(x[1],2) + \
3*x[0]*x[1]*x[2]*r_Z - 2*x[0]*Power(x[2],2)*r_Y - \
2*x[1]*Power(x[2],2)*r_X) - \
2*Power(a,2)*x[2]*Power(r,6)*(-3*a*x[1]*r_Y + x[0]*(x[1] + 3*a*r_X)) \
+ 4*a*Power(r,8)*((Power(x[0],2) - Power(x[1],2))*r_Z + \
x[2]*(x[1]*r_Y - x[0]*r_X)) + \
2*Power(a,2)*Power(r,7)*(3*x[0]*x[1]*r_Z + 2*x[2]*(x[0]*r_Y + \
x[1]*r_X)) + Power(r,9)*(-3*x[0]*x[1]*r_Z + 3*x[2]*(x[0]*r_Y + \
x[1]*r_X)) + Power(a,2)*Power(r,5)*(x[0]*x[1]*(Power(a,2) + \
Power(x[2],2))*r_Z + x[2]*(2*a*(Power(x[0],2) - Power(x[1],2)) + \
x[0]*(Power(a,2) - Power(x[2],2))*r_Y + x[1]*(Power(a,2) - \
Power(x[2],2))*r_X)))) - x[5]*x[6]*r*(2*G*M*Power(r,4)*Power(a*x[0] - \
x[1]*r,2)*(4*Power(r,2)*(a*x[0] - x[1]*r)*(Power(a,2)*Power(x[2],2) + \
Power(r,4))*r_Z + 2*x[1]*r*(Power(a,2) + \
Power(r,2))*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_Z - 3*(a*x[0] - \
x[1]*r)*(Power(a,2) + Power(r,2))*(Power(a,2)*Power(x[2],2) + \
Power(r,4))*r_Z + r*(a*x[0] - x[1]*r)*(Power(a,2) + \
Power(r,2))*(2*Power(a,2)*x[2] + 4*Power(r,3)*r_Z)) + \
2*G*M*Power(x[2],3)*r*Power(Power(a,2) + \
Power(r,2),4)*(Power(a,2)*Power(x[2],2) - 3*Power(r,4))*r_Y - \
Power(Power(a,2) + Power(r,2),2)*(Power(a,2)*Power(x[2],2) + \
2*G*M*Power(r,3) + Power(r,4))*(2*Power(r,3)*(a*x[0] - \
x[1]*r)*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_Z + \
x[1]*Power(r,2)*(Power(a,2) + Power(r,2))*(Power(a,2)*Power(x[2],2) + \
Power(r,4))*r_Z + 3*r*(-(a*x[0]) + x[1]*r)*(Power(a,2) + \
Power(r,2))*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_Z + \
Power(r,2)*(a*x[0] - x[1]*r)*(Power(a,2) + \
Power(r,2))*(2*Power(a,2)*x[2] + 4*Power(r,3)*r_Z) - \
4*x[2]*Power(r,4)*Power(Power(a,2) + Power(r,2),2)*r_Y + \
2*x[2]*Power(Power(a,2) + Power(r,2),2)*(Power(a,2)*Power(x[2],2) + \
Power(r,4))*r_Y) + 2*G*M*Power(r,3)*(a*x[1] + \
x[0]*r)*(2*Power(a,7)*Power(x[2],3)*(x[1]*r_Y + x[0]*r_X) - \
2*Power(a,5)*x[2]*Power(r,4)*(x[1]*r_Y + x[0]*r_X) + \
2*Power(a,5)*x[2]*Power(r,2)*(a*x[0]*x[1] + 2*(-Power(x[0],2) + \
Power(x[1],2))*x[2]*r_Z + x[1]*Power(x[2],2)*r_Y + \
x[0]*Power(x[2],2)*r_X) + \
2*Power(a,4)*x[2]*Power(r,3)*(a*Power(x[0],2) - a*Power(x[1],2) + \
2*a*Power(x[2],2) + 3*x[0]*x[1]*x[2]*r_Z + 2*x[0]*Power(x[2],2)*r_Y - \
2*x[1]*Power(x[2],2)*r_X) - \
2*Power(a,2)*x[2]*Power(r,6)*(3*a*x[1]*r_Y + x[0]*(x[1] + 3*a*r_X)) + \
4*a*Power(r,8)*((Power(x[0],2) - Power(x[1],2))*r_Z - x[2]*(x[1]*r_Y \
+ x[0]*r_X)) + Power(a,6)*Power(x[2],2)*r*(-3*x[0]*x[1]*r_Z + \
x[2]*(2*a + 3*x[0]*r_Y - 3*x[1]*r_X)) + \
2*Power(a,2)*Power(r,7)*(3*x[0]*x[1]*r_Z + 2*x[2]*(a - x[0]*r_Y + \
x[1]*r_X)) + Power(r,9)*(-3*x[0]*x[1]*r_Z + x[2]*(2*a - 3*x[0]*r_Y + \
3*x[1]*r_X)) + Power(a,2)*Power(r,5)*(x[0]*x[1]*(Power(a,2) + \
Power(x[2],2))*r_Z + x[2]*(2*a*(Power(a,2) + Power(x[0],2) - \
Power(x[1],2) + Power(x[2],2)) + x[0]*(-Power(a,2) + \
Power(x[2],2))*r_Y + x[1]*(Power(a,2) - Power(x[2],2))*r_X)))) + \
x[4]*x[6]*r*(2*G*M*Power(r,4)*Power(a*x[1] + \
x[0]*r,2)*(4*Power(r,2)*(a*x[1] + x[0]*r)*(Power(a,2)*Power(x[2],2) + \
Power(r,4))*r_Z - 2*x[0]*r*(Power(a,2) + \
Power(r,2))*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_Z - 3*(a*x[1] + \
x[0]*r)*(Power(a,2) + Power(r,2))*(Power(a,2)*Power(x[2],2) + \
Power(r,4))*r_Z + 2*r*(a*x[1] + x[0]*r)*(Power(a,2) + \
Power(r,2))*(Power(a,2)*x[2] + 2*Power(r,3)*r_Z)) - \
2*G*M*Power(x[2],3)*r*Power(Power(a,2) + \
Power(r,2),4)*(Power(a,2)*Power(x[2],2) - 3*Power(r,4))*r_X - \
Power(Power(a,2) + Power(r,2),2)*(Power(a,2)*Power(x[2],2) + \
2*G*M*Power(r,3) + Power(r,4))*(2*Power(r,3)*(a*x[1] + \
x[0]*r)*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_Z - \
x[0]*Power(r,2)*(Power(a,2) + Power(r,2))*(Power(a,2)*Power(x[2],2) + \
Power(r,4))*r_Z - 3*r*(a*x[1] + x[0]*r)*(Power(a,2) + \
Power(r,2))*(Power(a,2)*Power(x[2],2) + Power(r,4))*r_Z + \
2*Power(r,2)*(a*x[1] + x[0]*r)*(Power(a,2) + \
Power(r,2))*(Power(a,2)*x[2] + 2*Power(r,3)*r_Z) + \
4*x[2]*Power(r,4)*Power(Power(a,2) + Power(r,2),2)*r_X - \
2*x[2]*Power(Power(a,2) + Power(r,2),2)*(Power(a,2)*Power(x[2],2) + \
Power(r,4))*r_X) + 2*G*M*Power(r,3)*(a*x[0] - \
x[1]*r)*(-2*Power(a,7)*Power(x[2],3)*(x[1]*r_Y + x[0]*r_X) + \
2*Power(a,5)*x[2]*Power(r,4)*(x[1]*r_Y + x[0]*r_X) + \
2*Power(a,2)*x[2]*Power(r,6)*(-(x[0]*x[1]) + 3*a*x[1]*r_Y + \
3*a*x[0]*r_X) + 2*Power(a,5)*x[2]*Power(r,2)*(a*x[0]*x[1] + \
2*(-Power(x[0],2) + Power(x[1],2))*x[2]*r_Z - x[1]*Power(x[2],2)*r_Y \
- x[0]*Power(x[2],2)*r_X) + \
2*Power(a,4)*x[2]*Power(r,3)*(a*Power(x[0],2) - a*Power(x[1],2) - \
2*a*Power(x[2],2) + 3*x[0]*x[1]*x[2]*r_Z - 2*x[0]*Power(x[2],2)*r_Y + \
2*x[1]*Power(x[2],2)*r_X) + 4*a*Power(r,8)*((Power(x[0],2) - \
Power(x[1],2))*r_Z + x[2]*(x[1]*r_Y + x[0]*r_X)) - \
Power(a,6)*Power(x[2],2)*r*(3*x[0]*x[1]*r_Z + x[2]*(2*a + 3*x[0]*r_Y \
- 3*x[1]*r_X)) + 2*Power(a,2)*Power(r,7)*(3*x[0]*x[1]*r_Z - 2*x[2]*(a \
- x[0]*r_Y + x[1]*r_X)) - Power(r,9)*(3*x[0]*x[1]*r_Z + x[2]*(2*a - \
3*x[0]*r_Y + 3*x[1]*r_X)) + \
Power(a,2)*Power(r,5)*(x[0]*x[1]*(Power(a,2) + Power(x[2],2))*r_Z - \
x[2]*(2*a*(Power(a,2) - Power(x[0],2) + Power(x[1],2) + \
Power(x[2],2)) + x[0]*(-Power(a,2) + Power(x[2],2))*r_Y + \
x[1]*(Power(a,2) - Power(x[2],2))*r_X))))))/(Power(Power(a,2) + \
Power(r,2),4)*Power(Power(a,2)*Power(x[2],2) + Power(r,4),3));

        break;
    default:
        std::cerr << "invalid switch \n";
    }
   }



void next_ray_on_geodesic(ray& r, step& lambda, state_type& x,  runge_kutta4< state_type > rk4) {
    const auto r_old = r;
    double radius = (point3(x[0], x[1], x[2]) - star.location()).length();
    
    rk4.do_step(geodesicODE, x, lambda.t, lambda.adjusted_step_size(radius));
    update_ray(r_old, r, x);
    lambda.t += lambda.adjusted_step_size(radius);
}


color ray_color( ray& r, const hittable_list& world, step& lambda, const double background_radius, int bounce_depth, int step_depth, state_type& x, runge_kutta4< state_type > rk4) {

    hit_record rec;
    
    next_ray_on_geodesic(r, lambda, x, rk4);
    
    double radius = (point3(x[0], x[1], x[2]) - star.location()).length();
    
    //if we've exceeded the ray bounce limit, no more light is gathered.
    if (bounce_depth <= 0) {
        return color(0, 0, 0);
    }

    //if we exceed search length return pink for debugging (avoids infinte loop)
    if (step_depth <= 0) {
        return color(1, 0, 1);
    }

    //ray hits object
    if (world.hit(r, 0.00001, 1, rec)) {
        ray scattered;
        color attenuation;
        color emitted = rec.mat_ptr->emitted(rec.u, rec.v, rec.p);

        if (!rec.mat_ptr->scatter(r, rec, attenuation, scattered))
            if (emitted.e[3] < 1) {
                return emitted.e[3] * emitted + (1 - emitted.e[3]) * ray_color(r, world, lambda, background_radius, bounce_depth - 1, step_depth, x, rk4);
            }
            return emitted;
        
        
        
        if (rec.mat_ptr->scatter(r, rec, attenuation, scattered)) {
            if (attenuation == color(0, 0, 0) && attenuation.e[3] == 1.0)
                return attenuation;

            if (attenuation.e[3] < 1) {
                color passed_through = (1 - attenuation.e[3]) * ray_color(r, world, lambda, background_radius, bounce_depth - 1, step_depth, x, rk4);
                update_x(scattered, x);
                r.dir = rec.p - r.orig;
                return passed_through + attenuation.e[3] * attenuation * ray_color(r, world, lambda, background_radius, bounce_depth - 1, step_depth, x, rk4);
            }
            update_x(scattered, x);
            r.dir = rec.p - r.orig;
            return  6 * attenuation * ray_color(r, world, lambda, background_radius, bounce_depth - 1, step_depth, x, rk4);
        }
        return color(0, 0, 0);
    }

    //ray does not hit object and is out of bound
    if (radius > background_radius) {
        if (world.background_exits) {
            if (world.background_object->hit(ray(world.camera_position, unit_vector(r.dir)), 0.00001, infinity, rec)) {
                ray scattered;
                color attenuation;

                if (rec.mat_ptr->scatter(ray(world.camera_position, point3(x[4], x[5], x[6])), rec, attenuation, scattered)) {
                    return attenuation;
                }
            }
        }
        // no background object so give sky colours
        return color(0, 0, 0);
    }

    //ray does not hit object and is in bounds
    return ray_color(r, world, lambda, background_radius, bounce_depth, step_depth-1, x, rk4);
}


hittable_list earth(const point3 lookfrom) {
    hittable_list world(lookfrom);


    auto black_hole_texture = make_shared<solid_color>(color(0, 0, 0));
    auto black_hole_surface = make_shared<lambertian>(black_hole_texture);
    auto black_hole = make_shared<sphere>(star.location(), r_schwarzschild, black_hole_surface);
    world.add(black_hole);

    auto space_texture = make_shared<image_texture>("milkyway_2020_8k.jpg");
    auto space_surface = make_shared<lambertian>(space_texture);
    auto space = make_shared<sphere>(lookfrom, 500, space_surface);
    //world.add_background(space);

    const auto ring_factor = 8.2;
    auto ring_texture = make_shared<image_texture>("black_hole_ring.png");
    auto ring_material = make_shared<diffuse_light>(ring_texture);
    world.add(make_shared<xy_rect>(-1 * ring_factor* r_schwarzschild, 1 * ring_factor * r_schwarzschild, -1 * ring_factor * r_schwarzschild, 1 * ring_factor * r_schwarzschild, 0, ring_material));




    return hittable_list(world);
}


int main() {

    // Image
    const auto aspect_ratio = 9.0 / 9.0;
    const int image_width = 800;
    const int image_height = static_cast<int>(image_width / aspect_ratio);
    const int samples_per_pixel = 16;
    const int max_bounce_depth = 10;

    const int channel_num = 3;
    unsigned char* image_data = new unsigned char[image_width * image_height * channel_num];
    const std::string folder = "images/";
    const std::string ending = ".png";
    std::string filename = "image";
    std::string image_name = folder + filename + ending;


    // Geodesic

    const auto background_radius = 20 * r_schwarzschild;
    const auto max_geodesic_length = background_radius * 30;
    int geodesic_depth;
    switch (metric_case) 
    {
    case 0:                      // No curvature
        geodesic_depth = 10;
        break;
    default:                        
        geodesic_depth = 800;
        break;
    }

    step lambda;
    lambda.dt = static_cast<double>(max_geodesic_length / geodesic_depth);
    lambda.dt_min = lambda.dt / 4;
    lambda.dt_max = lambda.dt * 1;
    lambda.t = 0.0;
    lambda.use_adjusted_step_size = false;
    lambda.r_schwarzschild = r_schwarzschild;
    
    state_type x(8);
    runge_kutta4< state_type > rk4;

    // Camera

    point3 lookfrom(-12, -24, 5);
    point3 lookat(0.01, 0.01, 0.01);
    vec3 vup(0, 0, 1);
    double vof(90);

    camera cam(lookfrom, lookat, vup, vof, aspect_ratio);


    // World

    auto world = earth(lookfrom);
    
    world.add_frame(lookfrom, 0.99);

    //for (int i = 0; i < 8; i++) {
    //    lookfrom = point3(10 * sin(0.5 * pi * (i / 7.0)), 0, 10 * cos(0.5 * pi * (i / 7.0)));
    //    world.add_frame(lookfrom);
    //}

    //world.add_frame(lookfrom, 0.0);
    //world.add_frame(lookfrom ,0.2);
    //world.add_frame(lookfrom, 0.4);
    //world.add_frame(lookfrom, 0.6);
    //world.add_frame(lookfrom, 0.8);
    //world.add_frame(lookfrom, 0.9);
    //world.add_frame(lookfrom, 0.95);

    for (int frame = 0; frame < world.frame_num; frame++) {
        cam = camera(world.camera_position_frame(frame), lookat, vup, vof, aspect_ratio);
        world.camera_position = world.camera_position_frame(frame);
        a = world.rotation_frame(frame);
        
        filename = "image";
        std::cerr << "Frame (" << frame + 1 << "/" << world.frame_num << "): " << filename << "\n";
        for (int j = image_height - 1; j >= 0; --j) {
            std::cerr << "\rScanlines rendered: (" << image_height - j << '/' << image_height << ')' << std::flush;
            for (int i = 0; i < image_width; ++i) {
                color pixel_color(0, 0, 0);
                for (int s = 0; s < samples_per_pixel; s++) {
                    auto u = (i + random_double()) / (image_width - 1);
                    auto v = (j + random_double()) / (image_height - 1);
                    ray r = cam.get_ray(u, v);
                    x[3] = 0;
                    x[7] = -1;
                    update_x(r, x);
                    r.dir = r.orig;
                    pixel_color += ray_color(r, world, lambda, background_radius, max_bounce_depth, geodesic_depth, x, rk4);
                }
                write_color(image_data, pixel_color, samples_per_pixel, i, j, image_width, image_height, channel_num);
            }
        }
        stbi_write_png((folder + filename + ending).c_str(), image_width, image_height, channel_num, image_data, image_width * channel_num);
        std::cerr << "\n";
    }   
    // Render

    delete[] image_data;
    std::cerr << "\nDone.\n";
}


