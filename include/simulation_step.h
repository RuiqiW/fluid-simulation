#ifndef SIMULATION_STEP_H
#define SIMULATION_STEP_H

#include <Eigen/Dense>
#include <particle_property.h>
#include <collision_plane.h>
#include <vector>

void simulation_step(	
    Particles& particles, 
	const Eigen::MatrixXd &V_wall,
    const Eigen::MatrixXi &F_wall,
	const Eigen::MatrixXd &N_wall,
    double dt);

#endif
