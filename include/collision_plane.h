#ifndef COLLISION_PLANE_H
#define COLLISION_PLANE_H

#include <Eigen/Dense>
#include <particle_property.h>
#include <boundary_property.h>


void collision_plane(
    Particles& particles,
    Eigen::MatrixXd& pred_position,
    Wall_Plane &wall);


#endif
