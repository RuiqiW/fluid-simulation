#include <boundary_property.h>
#include <igl/AABB.h>
#include <igl/list_to_matrix.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/per_face_normals.h>
#include <igl/readOBJ.h>
#include <init_boundaries.h>
#include <init_position.h>
#include <particle_property.h>
#include <simulation_step.h>

#include <Eigen/Core>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
using namespace std;


// viewer window
igl::opengl::glfw::Viewer viewer;

// the boundary of box
// TODO: change color
void draw_boudary(const Eigen::MatrixXd& V_box) {

    // Edges of the bounding box
    viewer.data().add_points(V_box, Eigen::RowVector3d(1, 0, 0));
    Eigen::MatrixXi E_box(12, 2);
    E_box << 0, 1,
        1, 2,
        2, 3,
        3, 0,
        4, 5,
        5, 6,
        6, 7,
        7, 4,
        0, 4,
        1, 5,
        2, 6,
        7, 3;
    // Plot the edges of the bounding box
    for (unsigned i = 0; i < E_box.rows(); ++i) {
        viewer.data().add_edges(
            V_box.row(E_box(i, 0)),
            V_box.row(E_box(i, 1)),
            Eigen::RowVector3d(1, 0, 0));
    }
}

int main(int argc, char* argv[]) {
    if (argc <= 3 && argc >=2) {
        // choose whether it should be animation
        bool isAnimation;
        if (argv[1] == std::string("y")) {
            isAnimation = true;
        }else if (argv[1] == std::string("n")) {
            // choose whether it should be animation
            isAnimation = false;
        }
        else{
            std::cout << "Invalid argument to show the visualization , exiting.\n";
            exit(1);
        }
        const int xid = viewer.selected_data_index;
        viewer.append_mesh();

        // TODO: change the wording here
        std::cout << R"(
          N,n      Get to the next timestep for fluid simulation. 
          R,r      Reset all particles' positions and velocities.
        )";

        //about position and velocity
        //starting corner for the fluid particles;
        Particles particles;
        Eigen::Vector3d corner(-1., -0.2, -1.);
        Eigen::Vector3i num_points(20, 20, 20);
        double step_size = 0.054;
        const Eigen::RowVector3d particle_color(0.1, 0.9, 0.9);

        init_position(particles, corner, num_points, step_size);

        /* Add Walls */
        std::cout << "Add Walls" << std::endl;
        Eigen::Vector3d dm, dM;
        dm << 0.2, 1.0, 1.0;
        dM << 1.5, 0.5, 1.0;
        Eigen::Vector3d m = particles.position.colwise().minCoeff().transpose() - dm;
        Eigen::Vector3d M = particles.position.colwise().maxCoeff().transpose() + dM;

        Rabbit_Mesh rabbit;
        if (argc == 3) {
            if (argv[2] == std::string("rabbit")) {
                /* Add Bunny */
                init_rabbit(rabbit, "../data/coarser_bunny.obj", true);
                // init_rabbit(rabbit);
                viewer.data().set_mesh(rabbit.V_obj, rabbit.F_obj);
                m(1) = rabbit.V_obj.col(1).minCoeff();

            } else {
                std::cout << "Invalid argument to show the visualization , exiting.\n";
                exit(1);
            }

        }

        else if (argc == 2) {
            init_rabbit(rabbit, "", false);
        }

        Wall_Plane wall;
        init_wall(wall, m, M, true);

        igl::per_face_normals(wall.V_wall, wall.F_wall, wall.N_wall);

        const auto& reset = [&]() {
            init_position(particles, corner, num_points, step_size);
            viewer.data_list[xid].set_points(particles.position, (1. - (1. - particle_color.array()) * .9));
        };

        const auto& move = [&]() {
            // TODO: fill the algorithm
            simulation_step(particles, wall, rabbit, double(0.016));
            viewer.data_list[xid].set_points(particles.position, (1. - (1. - particle_color.array()) * .9));
        };

        if (isAnimation == false) {
            viewer.callback_key_pressed =
                [&](igl::opengl::glfw::Viewer&, unsigned char key, int) -> bool {
                switch (key) {
                    case 'N':
                    case 'n':

                        move();
                        break;
                    case 'R':
                    case 'r':
                        reset();
                        break;
                    default:
                        return false;
                }
                return true;
            };
        }

        std::cout << "draw_boudary" << std::endl;
        // draw_boudary(wall.V_wall);
        std::cout << "viewer" << std::endl;

        viewer.data_list[xid].set_colors(particle_color);
        viewer.data_list[xid].set_points(particles.position, (1. - (1. - particle_color.array()) * .9));
        viewer.data_list[xid].point_size = 6.0;

        if (isAnimation) {
            viewer.callback_pre_draw = [&](igl::opengl::glfw::Viewer&) -> bool {
                // Create orbiting animation
                simulation_step(particles, wall, rabbit, double(0.016));
                viewer.data_list[xid].set_points(particles.position, (1. - (1. - particle_color.array()) * .9));
                return false;
            };
            
        }
        viewer.launch();
    }
    return EXIT_SUCCESS;
}