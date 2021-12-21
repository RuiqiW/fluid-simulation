#ifndef BOUNDARY_PROPERTY_H
#define BOUNDARY_PROPERTY_H

#include <Eigen/Dense>
#include <igl/AABB.h>

using namespace std;

struct Wall_Plane {
    bool isUsed;
    Eigen::MatrixXd V_wall; // Corners of the bounding box
    Eigen::MatrixXd N_wall; // normal of the each face that wall consists of 
    Eigen::MatrixXi F_wall; 

};

struct Rabbit_Mesh {
    bool isUsed;
    Eigen::MatrixXd V_obj; // vertex of rabbit
    Eigen::MatrixXi F_obj; // face of rabbit
    igl::AABB<Eigen::MatrixXd, 3> tree;

};

#endif