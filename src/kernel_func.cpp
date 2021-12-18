#include "kernel_func.h"
#include "params.h"
#include <math.h>

double POLY6_COEF = 315.0 / (64.0 * M_PI * std::pow(RADIUS, 9));

double GRAD_POLY_COEF = -945.0 / (32.0 * M_PI * std::pow(RADIUS, 9));

double poly6(double r){
    return std::pow(RADIUS * RADIUS - r * r, 3) * POLY6_COEF;
}

void dPoly6(Eigen::Vector3d &dW, Eigen::Ref<const Eigen::Vector3d> d, double r){
    dW = GRAD_POLY_COEF * std::pow(RADIUS  * RADIUS - r * r, 2) * d;
}


double SPIKY_COEF = 15.0 / (M_PI * std::pow(RADIUS, 6));

double GRAD_SPIKY_COEF = - 45.0 / (M_PI * std::pow(RADIUS, 6));


double spiky(double r){
    return SPIKY_COEF * std::pow(RADIUS - r, 3);
}

void dSpiky(Eigen::Vector3d &dW, Eigen::Ref<const Eigen::Vector3d> d, double r){
    if(r < 1e-7){
        dW = Eigen::Vector3d::Zero();
    }else{
        dW = GRAD_SPIKY_COEF * std::pow(RADIUS - r, 2) / r * d;
    }
}

