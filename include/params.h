#ifndef PARAMS_H
#define PARAMS_H
#include <Eigen/Dense>

#define RADIUS 0.1
#define ITERS 4
#define RHO_REST 6378.0
#define MASS 1.0
#define CFM_EPS 600.0

// artificial pressure
#define P_RADIUS 0.03
#define P_POWER 4.0
#define P_CONST 0.0001

#define GRAVITY Eigen::Vector3d(0.0, -9.8, 0.0)

//for XSPH viscosity
#define C_VISCOSITY 0.0001
#define EPS_VORTICITY 0.000000001
#endif