#include "kernel_func.h"
#include "params.h"
#include <math.h>

double POLY6_COEF = 315.0 / (64.0 * M_PI * std::pow(RADIUS, 9));

double GRAD_POLY_COEF = -945.0 / (32.0 * M_PI * std::pow(RADIUS, 9));

double poly6(double r){
    return std::pow(RADIUS * RADIUS - r * r, 3) * POLY6_COEF;
}

void dPoly6(Eigen::Vector3d &dC, Eigen::Ref<const Eigen::Vector3d> d, double r){
    dC = GRAD_POLY_COEF * std::pow(RADIUS  * RADIUS - r * r, 2) * d;
}

