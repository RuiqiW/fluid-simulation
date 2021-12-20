#ifndef SIMULATION_STEP_H
#define SIMULATION_STEP_H

#include <Eigen/Dense>
#include <particle_property.h>
#include <collision_plane.h>
#include <collision_mesh.h>
#include <vector>
#include <igl/AABB.h>
#include <boundary_property.h>

void simulation_step(
    Particles &particles,
    Wall_Plane &wall,
    Rabbit_Mesh &rabbit,
    double dt);

#endif
