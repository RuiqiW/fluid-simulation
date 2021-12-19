#include <kernel_func.h>
#include <params.h>
#include <simulation_step.h>
#include <find_neighbors.h>

#include <iostream>
#include <vector>


void simulation_step(
    Particles &particles,
    Wall_Plane &wall,
    Rabbit_Mesh &rabbit,
    double dt) {

    int N = particles.position.rows();
    


    Eigen::MatrixXd pred_position;
    pred_position.resize(N, 3);

    // apply forces and prediction positions
    particles.velocity += dt * particles.acceleration;
    pred_position = particles.position + dt * particles.velocity;


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

    for (int i = 0; i < N; i++) {
        Eigen::Vector3d pos = pred_position.row(i);
        int x = (pos(0) - min_x) / RADIUS;
        int y = (pos(1) - min_y) / RADIUS;
        int z = (pos(2) - min_z) / RADIUS;

        Grid[x][y][z].push_back(i);
    }

    // TODO: find neighbours
    // find_neighbors(particles, pred_position);

    Eigen::MatrixXd pos_correction;
    pos_correction.resize(N, 3);

    Eigen::VectorXd lambdas;
    lambdas.resize(N);

    for (int iter = 0; iter < ITERS; iter++) {

        Eigen::VectorXd dC2;
        dC2.resize(N);
        dC2.setZero();

        // loop over grid
        for (int x = 0; x < dim_x; x++) {
            for (int y = 0; y < dim_y; y++) {
                for (int z = 0; z < dim_z; z++) {
                    for (int index : Grid[x][y][z]) {
                        double density = 0;  // Rho

                        Eigen::Vector3d gradI;  // dW / dpi
                        gradI.setZero();

                        double sum = 0;  // sum |dW / dpj|^2

                        for (int a = -1; a < 2; a++) {
                            if (x + a >= 0 && x + a < dim_x) {
                                for (int b = -1; b < 2; b++) {
                                    if (y + b >= 0 && y + b < dim_y) {
                                        for (int c = -1; c < 2; c++) {
                                            if (z + c >= 0 && z + c < dim_z) {
                                                for (int p : Grid[x + a][y + b][z + c]) {
                                                    if (p == index) continue;
                                                    Eigen::Vector3d d = pred_position.row(index) - pred_position.row(p);
                                                    double r = d.norm();
                                                    if (r < RADIUS) {
                                                        density += poly6(r);
                                                        
                                                        Eigen::Vector3d gradJ;
                                                        dSpiky(gradJ, d, r);
                                                        gradJ /= RHO_REST;
                                                        gradI += gradJ;
                                                        sum += gradJ.squaredNorm();
                                                        
                                                    }
                                                }
                                            }
                                        }  //end c
                                    }
                                }  //end b
                            }
                        }  // end a

                        particles.density(index) = density + poly6(0.0);
                        dC2(index) = sum + gradI.squaredNorm();
                    }  // end iteration over neighbors of a particle
                }  // end z
            }  // end y
        }  // end x


        /* calculate lambda */
        for (int i = 0; i < N; i++) {
            double Ci = std::max(particles.density(i) * 1.0 / RHO_REST - 1.0, 0.0);
            lambdas(i) = -1.0 * Ci / (dC2(i) + CFM_EPS);
        }


        double W_delta_q = poly6(P_RADIUS);

        // calculate position correction
        for (int x = 0; x < dim_x; x++) {
            for (int y = 0; y < dim_y; y++) {
                for (int z = 0; z < dim_z; z++) {
                    for (int index : Grid[x][y][z]) {
                        Eigen::Vector3d diff_pos;
                        diff_pos.setZero();

                        for (int a = -1; a < 2; a++) {
                            if (x + a >= 0 && x + a < dim_x) {
                                for (int b = -1; b < 2; b++) {
                                    if (y + b >= 0 && y + b < dim_y) {
                                        for (int c = -1; c < 2; c++) {
                                            if (z + c >= 0 && z + c < dim_z) {
                                                for (int p : Grid[x + a][y + b][z + c]) {
                                                    if (p == index) continue;
                                                    Eigen::Vector3d d = pred_position.row(index) - pred_position.row(p);
                                                    double r = d.norm();
                                                    if (r < RADIUS) {
                                                        double lmbd = lambdas(index) + lambdas(p);
                                                        Eigen::Vector3d gradW;
                                                        dSpiky(gradW, d, r);

                                                        double W = poly6(r);
                                                        double scorr = -P_CONST * std::pow(W / W_delta_q, P_POWER);
                                                        diff_pos += (lmbd + scorr) * gradW;
                                                    }
                                                }
                                            }
                                        }  //end c
                                    }
                                }  //end b
                            }
                        }  // end a

                        diff_pos /=  RHO_REST;
                        pos_correction.row(index) = diff_pos;
                    }  // end iteration over neighbors of a particle
                }  // end z
            }  // end y
        }  // end x

        // update predicted position
        pred_position += pos_correction;              


        // TODO: collision detection and response
        if (rabbit.isUsed){
            collision_mesh(particles, pred_position, rabbit);
        }
        if (wall.isUsed){
            collision_plane(particles, pred_position, wall);
        }
        
    }

    // update velocity
    particles.velocity = (pred_position - particles.position) / dt;    
    
    // update position
    particles.position = pred_position;

    
    /* apply vorticity confinement and XSPH viscosity */
    Eigen::MatrixXd vorticity;
    vorticity.resize(N, 3);

    Eigen::MatrixXd velocity_change;
    velocity_change.resize(N, 3);

    for (int x = 0; x < dim_x; x++) {
        for (int y = 0; y < dim_y; y++) {
            for (int z = 0; z < dim_z; z++) {
                for (int index : Grid[x][y][z]) {
                    
                    Eigen::Vector3d delta_v;
                    delta_v.setZero();

                    Eigen::Vector3d omega;
                    omega.setZero();

                    for (int a = -1; a < 2; a++) {
                        if (x + a >= 0 && x + a < dim_x) {
                            for (int b = -1; b < 2; b++) {
                                if (y + b >= 0 && y + b < dim_y) {
                                    for (int c = -1; c < 2; c++) {
                                        if (z + c >= 0 && z + c < dim_z) {
                                            for (int p : Grid[x + a][y + b][z + c]) {
                                                if (p == index) continue;
                                                Eigen::Vector3d d = pred_position.row(index) - pred_position.row(p);
                                                double r = d.norm();
                                                if (r < RADIUS) {
                                                    Eigen::Vector3d vij;
                                                    vij = particles.velocity.row(p) - particles.velocity.row(index);

                                                    //update change amount in velocity
                                                    double W = poly6(r);
                                                    delta_v +=  W * vij;

                                                    Eigen::Vector3d gradW;
                                                    dSpiky(gradW, d, r);
                                                    omega += vij.cross(gradW);
                                                }
                                            }
                                        }
                                    }  //end c
                                }
                            }  //end b
                        }
                    }  // end a

                    velocity_change.row(index) = C_VISCOSITY * delta_v;
                    vorticity.row(index) = omega;

                }  // end iteration over neighbors of a particle
            }  // end z
        }   // end y
    }    // end x


    particles.velocity += velocity_change;

    // std::cout << particles.density.minCoeff() << " " << particles.density.maxCoeff() << "\n";

    for (int x = 0; x < dim_x; x++) {
        for (int y = 0; y < dim_y; y++) {
            for (int z = 0; z < dim_z; z++) {
                for (int index : Grid[x][y][z]) {
                    Eigen::Vector3d omega = vorticity.row(index);
                    double omega_norm = omega.norm();

                    Eigen::Vector3d eta;
                    eta.setZero();
                    
                    for (int a = -1; a < 2; a++) {
                        if (x + a >= 0 && x + a < dim_x) {
                            for (int b = -1; b < 2; b++) {
                                if (y + b >= 0 && y + b < dim_y) {
                                    for (int c = -1; c < 2; c++) {
                                        if (z + c >= 0 && z + c < dim_z) {
                                            for (int p : Grid[x + a][y + b][z + c]) {
                                                if (p == index) continue;
                                                Eigen::Vector3d d = pred_position.row(index) - pred_position.row(p);
                                                double r = d.norm();
                                                if (r < RADIUS) {
                                                    Eigen::Vector3d gradW;
                                                    dSpiky(gradW, d, r);
                                                    eta += gradW * omega_norm;
                                                    
                                                }
                                            }
                                        }
                                    }  //end c
                                }
                            }  //end b
                        }
                    }  // end a

                    particles.acceleration.row(index) = GRAVITY;
                    Eigen::Vector3d f_vorticity;

                    if(eta.norm() != 0){
                        f_vorticity = EPS_VORTICITY * (eta.normalized()).cross(omega);
                        particles.acceleration.row(index) += f_vorticity / particles.density(index);
                    }

                }  // end iteration over neighbors of a particle
            }  // end z
        }   // end y
    }    // end x

}