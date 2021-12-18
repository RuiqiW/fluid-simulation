#include <particle_property.h>

#include <Eigen/Dense>
#include <iostream>
#include <params.h>
using namespace std;

void find_neighbors(Particles &particles, Eigen::Ref<const Eigen::MatrixXd> pred_position) {
    int N = particles.position.rows();

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

    vector<int> Grid[dim_x][dim_y][dim_z];

    for (int i = 0; i < N; i++) {
        Eigen::Vector3d pos = pred_position.row(i);
        int x = (pos(0) - min_x) / RADIUS;
        int y = (pos(1) - min_y) / RADIUS;
        int z = (pos(2) - min_z) / RADIUS;

        Grid[x][y][z].push_back(i);
    }

    // loop over grid
    for (int x = 0; x < dim_x; x++) {
        for (int y = 0; y < dim_y; y++) {
            for (int z = 0; z < dim_z; z++) {
                for (int index : Grid[x][y][z]) {
                    for (int a = -1; a < 2; a++) {
                        if (x + a >= 0 && x + a < dim_x) {
                            for (int b = -1; b < 2; b++) {
                                if (y + b >= 0 && y + b < dim_y) {
                                    for (int c = -1; c < 2; c++) {
                                        if (z + c >= 0 && z + c < dim_z) {
                                            for (int p : Grid[x + a][y + b][z + c]) {
                                                if (p == index) continue;
                                                particles.neighbours[index].push_back(p);
                                            }
                                        }
                                    }  //end c
                                }
                            }  //end b
                        }
                    }  // end a

                }  // end iteration over neighbors of a particle
            }      // end z
        }          // end y
    }              // end x
}
