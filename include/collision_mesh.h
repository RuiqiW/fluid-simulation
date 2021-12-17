#ifndef COLLISION_MESH_H
#define COLLISION_MESH_H

#include <Eigen/Dense>
#include <particle_property.h>
#include <igl/AABB.h>


void collision_mesh(
    Particles& particles,
    Eigen::MatrixXd& pred_position,
    const Eigen::MatrixXd &V_obj,
    const Eigen::MatrixXi &F_obj,
    igl::AABB<Eigen::MatrixXd, 3> &tree
);

#endif
