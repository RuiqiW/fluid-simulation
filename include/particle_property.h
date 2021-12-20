#ifndef PARTICLE_PROPERTY_H
#define PARTICLE_PROPERTY_H

#include <Eigen/Dense>
#include <iostream>
#include <vector>
using namespace std;

struct Particles {
    Eigen::MatrixXd position;
    Eigen::MatrixXd velocity;
    Eigen::MatrixXd acceleration;
    // TODO: might need to change the property it has
    Eigen::VectorXd density;

    vector<vector<int> > neighbours;

};

#endif