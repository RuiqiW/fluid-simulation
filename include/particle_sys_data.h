#ifndef PARTICLES_DATA_H
#define PARTICLES_DATA_H

#include <Eigen/Dense>

struct Particles {
	Eigen::MatrixXd position;
	Eigen::MatrixXd velocity;
	Eigen::MatrixXd acceleration;
	Eigen::VectorXd pressure;
	Eigen::VectorXd density;
};

#endif