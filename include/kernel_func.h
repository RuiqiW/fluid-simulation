#ifndef KERNEL_FUNC_H
#define KERNEL_FUNC_H

#include <Eigen/Dense>
#include <vector>

double poly6(double r);

void dPoly6(Eigen::Vector3d &dC, Eigen::Ref<const Eigen::Vector3d> d, double r);


double spiky(double r);

void dSpiky(Eigen::Vector3d &dW, Eigen::Ref<const Eigen::Vector3d> d, double r);

#endif
