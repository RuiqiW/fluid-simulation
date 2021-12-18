#ifndef SIMULATION_STEP_H
#define SIMULATION_STEP_H

#include <Eigen/Dense>
#include <particle_property.h>
#include <collision_plane.h>
#include <collision_mesh.h>
#include <vector>
#include <igl/AABB.h>
#include <boundary_property.h>

// void simulation_step(	
//     Particles& particles, 
// 	const Eigen::MatrixXd &V_wall,
//     const Eigen::MatrixXi &F_wall,
// 	const Eigen::MatrixXd &N_wall,
//     const Eigen::MatrixXd &V_obj,
//     const Eigen::MatrixXi &F_obj,
//     igl::AABB<Eigen::MatrixXd, 3> &tree,
//     double dt);

void simulation_step(
    Particles &particles,
    Wall_Plane &wall,
    Rabbit_Mesh &rabbit,
    double dt);

#endif
