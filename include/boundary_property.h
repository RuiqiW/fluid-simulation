#ifndef BOUNDARY_PROPERTY_H
#define BOUNDARY_PROPERTY_H

#include <Eigen/Dense>
#include <igl/AABB.h>

using namespace std;

struct Wall_Plane {
    bool isUsed;
    Eigen::MatrixXd V_wall; // Corners of the bounding box
    Eigen::MatrixXd N_wall;
    Eigen::MatrixXi F_wall;

};

struct Rabbit_Mesh {
    bool isUsed;
    Eigen::MatrixXd V_obj;
    Eigen::MatrixXi F_obj; 
    igl::AABB<Eigen::MatrixXd, 3> tree;

};

#endif