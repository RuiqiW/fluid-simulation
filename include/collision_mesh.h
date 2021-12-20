#ifndef COLLISION_MESH_H
#define COLLISION_MESH_H

#include <Eigen/Dense>
#include <particle_property.h>
#include <igl/AABB.h>
#include <boundary_property.h>


void collision_mesh(
    Particles& particles,
    Eigen::MatrixXd& pred_position,
    Rabbit_Mesh &rabbit
);





#endif
