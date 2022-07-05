#ifndef VEC3_H
#define VEC3_H

#include <cmath>
#include <iostream>

using std::sqrt;

class vec3 {
public:
    vec3() : e{ 0,0,0 } {}
    vec3(double e0, double e1, double e2) : e{ e0, e1, e2 } {}

    double x() const { return e[0]; }
    double y() const { return e[1]; }
    double z() const { return e[2]; }

    vec3 operator-() const { return vec3(-e[0], -e[1], -e[2]); }
    double operator[](int i) const { return e[i]; }
    double& operator[](int i) { return e[i]; }

    vec3& operator+=(const vec3& v) {
        e[0] += v.e[0];
        e[1] += v.e[1];
        e[2] += v.e[2];
        return *this;
    }

    vec3& operator*=(const double t) {
        e[0] *= t;
        e[1] *= t;
        e[2] *= t;
        return *this;
    }

    vec3& operator/=(const double t) {
        return *this *= 1 / t;
    }

    double length() const {
        return sqrt(length_squared());
    }

    double length_squared() const {
        return e[0] * e[0] + e[1] * e[1] + e[2] * e[2];
    }

    inline static vec3 random() {
        return vec3(random_double(), random_double(), random_double());
    }
    
    inline static vec3 random(double min, double max) {
        return vec3(random_double(min, max), random_double(min, max), random_double(min, max));
    }
    
    bool near_zero() const {
        // Return true if the vector is close to zero in all dimensions.
        const auto s = 1e-8;
        return (fabs(e[0]) < s) && (fabs(e[1]) < s) && (fabs(e[2]) < s);
    }

public:
    double e[3];
};

// vec3 Utility Functions
inline bool operator==(const vec3& u, const vec3& v) {
    return (u.e[0] == v.e[0]) && (u.e[1] == v.e[1]) && (u.e[2] == v.e[2]);
}

inline std::ostream& operator<<(std::ostream& out, const vec3& v) {
    return out << v.e[0] << ' ' << v.e[1] << ' ' << v.e[2];
}

inline vec3 operator+(const vec3& u, const vec3& v) {
    return vec3(u.e[0] + v.e[0], u.e[1] + v.e[1], u.e[2] + v.e[2]);
}

inline vec3 operator-(const vec3& u, const vec3& v) {
    return vec3(u.e[0] - v.e[0], u.e[1] - v.e[1], u.e[2] - v.e[2]);
}

inline vec3 operator*(const vec3& u, const vec3& v) {
    return vec3(u.e[0] * v.e[0], u.e[1] * v.e[1], u.e[2] * v.e[2]);
}

inline vec3 operator*(double t, const vec3& v) {
    return vec3(t * v.e[0], t * v.e[1], t * v.e[2]);
}

inline vec3 operator*(const vec3& v, double t) {
    return t * v;
}

inline vec3 operator/(vec3 v, double t) {
    return (1 / t) * v;
}

inline double dot(const vec3& u, const vec3& v) {
    return u.e[0] * v.e[0]
        + u.e[1] * v.e[1]
        + u.e[2] * v.e[2];
}

inline vec3 cross(const vec3& u, const vec3& v) {
    return vec3(u.e[1] * v.e[2] - u.e[2] * v.e[1],
        u.e[2] * v.e[0] - u.e[0] * v.e[2],
        u.e[0] * v.e[1] - u.e[1] * v.e[0]);
}

inline vec3 unit_vector(vec3 v) {
    return v / v.length();
}



vec3 random_in_unit_sphere() {
    while (true) {
        auto p = vec3::random(-1, 1);
        if (p.length_squared() >= 1) {
            continue;
        }
        return p;
    }
}

vec3 random_unit_vector() {
    return unit_vector(random_in_unit_sphere());
}


class color {
public:
    color() : e{ 0,0,0 } {}
    color(double e0, double e1, double e2) : e{e0, e1, e2, 1.0} {}
    color(double e0, double e1, double e2, double e3) : e{ e0, e1, e2, e3 } {}

    double x() const { return e[0]; }
    double y() const { return e[1]; }
    double z() const { return e[2]; }
    double a() const { return e[3]; }

    color operator-() const { return color(-e[0], -e[1], -e[2]); }
    double operator[](int i) const { return e[i]; }
    double& operator[](int i) { return e[i]; }

    color& operator+=(const color& v) {
        e[0] += v.e[0];
        e[1] += v.e[1];
        e[2] += v.e[2];
        return *this;
    }

    color& operator*=(const double t) {
        e[0] *= t;
        e[1] *= t;
        e[2] *= t;
        return *this;
    }

    color& operator/=(const double t) {
        return *this *= 1 / t;
    }

    double length() const {
        return sqrt(length_squared());
    }

    double length_squared() const {
        return e[0] * e[0] + e[1] * e[1] + e[2] * e[2];
    }

    inline static color random() {
        return color(random_double(), random_double(), random_double());
    }

    inline static color random(double min, double max) {
        return color(random_double(min, max), random_double(min, max), random_double(min, max));
    }

    bool near_zero() const {
        // Return true if the vector is close to zero in all dimensions.
        const auto s = 1e-8;
        return (fabs(e[0]) < s) && (fabs(e[1]) < s) && (fabs(e[2]) < s);
    }

public:
    double e[4];
};

// color Utility Functions
inline bool operator==(const color& u, const color& v) {
    return (u.e[0] == v.e[0]) && (u.e[1] == v.e[1]) && (u.e[2] == v.e[2]);
}

inline std::ostream& operator<<(std::ostream& out, const color& v) {
    return out << v.e[0] << ' ' << v.e[1] << ' ' << v.e[2];
}

inline color operator+(const color& u, const color& v) {
    return color(u.e[0] + v.e[0], u.e[1] + v.e[1], u.e[2] + v.e[2]);
}

inline color operator-(const color& u, const color& v) {
    return color(u.e[0] - v.e[0], u.e[1] - v.e[1], u.e[2] - v.e[2]);
}

inline color operator*(const color& u, const color& v) {
    return color(u.e[0] * v.e[0], u.e[1] * v.e[1], u.e[2] * v.e[2]);
}

inline color operator*(double t, const color& v) {
    return color(t * v.e[0], t * v.e[1], t * v.e[2]);
}

inline color operator*(const color& v, double t) {
    return t * v;
}

inline color operator/(color v, double t) {
    return (1 / t) * v;
}


// Type aliases for vec3
using point3 = vec3;   // 3D point

#endif
