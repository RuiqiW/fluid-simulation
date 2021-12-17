#ifndef COLLISION_PLANE_H
#define COLLISION_PLANE_H

#include <Eigen/Dense>
#include <particle_property.h>


void collision_plane(
    Particles& particles,
    Eigen::MatrixXd& pred_position,
    const Eigen::MatrixXd &V_wall,
    const Eigen::MatrixXi &F_wall,
	const Eigen::MatrixXd &N_wall);

#endif
