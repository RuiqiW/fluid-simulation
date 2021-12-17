#include <vector>
#include <iostream>
#include <params.h>
#include <simulation_step.h>
#include <kernel_func.h>


void simulation_step(    
    Particles& particles, 	
	const Eigen::MatrixXd &V_wall,
    const Eigen::MatrixXi &F_wall,
	const Eigen::MatrixXd &N_wall,
    const Eigen::MatrixXd &V_obj,
    const Eigen::MatrixXi &F_obj,
    igl::AABB<Eigen::MatrixXd, 3> &tree,
    double dt) {

    int N = particles.position.rows();

    Eigen::MatrixXd pred_position;
    pred_position.resize(N, 3);

    // apply forces and prediction positions
    for(int i=0; i < N; i++){
        particles.velocity.row(i) += dt * particles.acceleration.row(i);
		pred_position.row(i) = particles.position.row(i) + dt * particles.velocity.row(i);
    }


    // binning particles based on predicted position

    double min_x = pred_position.col(0).minCoeff();
	double min_y = pred_position.col(1).minCoeff();
	double min_z = pred_position.col(2).minCoeff();

	double max_x = pred_position.col(0).maxCoeff();
	double max_y = pred_position.col(1).maxCoeff();
	double max_z = pred_position.col(2).maxCoeff();

    int dim_x = (max_x - min_x) / RADIUS + 1;
    int dim_y = (max_y - min_y) / RADIUS + 1;
    int dim_z = (max_z - min_z) / RADIUS + 1;

    std::vector<int> Grid[dim_x][dim_y][dim_z];

    for(int i=0; i < N; i++) {
        Eigen::Vector3d pos = pred_position.row(i);
        int x = (pos(0) - min_x) / RADIUS;
        int y = (pos(1) - min_y) / RADIUS;
        int z = (pos(2) - min_z) / RADIUS;

        Grid[x][y][z].push_back(i);
    }



    Eigen::MatrixXd pos_correction;
    pos_correction.resize(N, 3);

    Eigen::VectorXd lambdas;
    lambdas.resize(N);


    for (int iter = 0; iter < ITERS; iter++){
        Eigen::VectorXd densities;
        densities.resize(N);
        densities.setZero();

        // sum of dW^2 for neighboring particles
        Eigen::VectorXd dWj2;
        dWj2.resize(N);
        dWj2.setZero();

        // sum of dW for particle i
        Eigen::MatrixXd dW;
        dW.resize(N, 3);
        dW.setZero();

        // loop over grid
        for (int x=0; x < dim_x; x++){
            for(int y=0; y < dim_y; y++){
                for(int z=0; z < dim_z; z++){
                    
                    for(int index: Grid[x][y][z]){

                        double density = 0; // Rho

                        Eigen::Vector3d grad; // dW / dpi
                        grad.setZero();

                        double sqr_norm_sum = 0; // sum |dW / dpj|^2


                        for(int a=-1; a<2; a++){
                            if(x+a >= 0 && x+a < dim_x){

                                for(int b=-1; b<2; b++){
                                    if(y+b >=0 && y+b < dim_y){
                
                                        for(int c=-1; c<2; c++){
                                            if(z+c >=0 && z+c < dim_z){

                                                for(int p : Grid[x+a][y+b][z+c]){
                                                    if (p == index) continue;
                                                    Eigen::Vector3d d = pred_position.row(index) - pred_position.row(p);
                                                    double r= d.norm();
                                                    if(r < RADIUS){
                                                        density += poly6(r);
                                                        Eigen::Vector3d tmp;
                                                        dPoly6(tmp, d, r);
                                                        grad += tmp;
                                                        sqr_norm_sum += tmp.squaredNorm();
                                                    }
                                                }
                                            }
                                        } //end c
                                    }
                                } //end b
                            }
                        } // end a

                        densities(index) = density;
                        dWj2(index) = sqr_norm_sum;
                        dW.row(index) = grad;
                    } // end iteration over neighbors of a particle
                } // end z
            } // end y
        } // end x

        
        // calculate lambda
        for(int i=0; i < N; i++){
            double Ci = densities(i) * 1.0 / RHO_REST - 1.0;

            double dC2 = (dWj2(i) + dW.row(i).squaredNorm()) / std::pow(RHO_REST, 2.0);
            // add relaxation parameter
            dC2 += CFM_EPS;

            lambdas(i) = - Ci / dC2;
        }

        double W_delta_q = poly6(P_RADIUS);

        // calculate position correction
        for (int x=0; x < dim_x; x++){
            for(int y=0; y < dim_y; y++){
                for(int z=0; z < dim_z; z++){
                    
                    for(int index: Grid[x][y][z]){

                        Eigen::Vector3d diff_pos;
                        diff_pos.setZero();

                        for(int a=-1; a<2; a++){
                            if(x+a >= 0 && x+a < dim_x){

                                for(int b=-1; b<2; b++){
                                    if(y+b >=0 && y+b < dim_y){
                
                                        for(int c=-1; c<2; c++){
                                            if(z+c >=0 && z+c < dim_z){

                                                for(int p : Grid[x+a][y+b][z+c]){
                                                    if (p == index) continue;
                                                    Eigen::Vector3d d = pred_position.row(index) - pred_position.row(p);
                                                    double r= d.norm();
                                                    if(r < RADIUS){
                                                        double lmbd = lambdas(index) + lambdas(p);
                                                        Eigen::Vector3d gradW;
                                                        dPoly6(gradW, d, r);

                                                        double W = poly6(r);
                                                        double scorr = - P_CONST * std::pow(W / W_delta_q, P_POWER);
                                                        diff_pos += (lmbd + scorr) * gradW;
                                                    }
                                                }
                                            }
                                        } //end c
                                    }
                                } //end b
                            }
                        } // end a

                        diff_pos = (1.0 / RHO_REST) * diff_pos;
                        pos_correction.row(index) = diff_pos;
                    } // end iteration over neighbors of a particle
                } // end z
            } // end y
        } // end x   

        // update predicted position
        pred_position += pos_correction;

        // TODO: collision detection and response
        collision_mesh(particles, pred_position, V_obj, F_obj, tree);
        collision_plane( particles, pred_position, V_wall, F_wall, N_wall);
    }





    // update velocity
    particles.velocity = (1.0 / dt) * (pred_position - particles.position);

    for(int i=0; i < N; i++){
        // TODO: apply voriticity confinement and XSPH viscosity
    }

    // update position
    particles.position = pred_position;
}