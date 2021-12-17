#include <igl/list_to_matrix.h>
#include <igl/opengl/glfw/Viewer.h>
#include <Eigen/Core>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cstdlib>

#include <igl/opengl/glfw/Viewer.h>
#include <init_position.h>
#include <particle_property.h>
#include <simulation_step.h>


int main(int argc, char* argv[])
{
    igl::opengl::glfw::Viewer viewer;
    const int xid = viewer.selected_data_index;
    viewer.append_mesh();

    // TODO: change the wording here
    std::cout << R"(
          A,a      Move one timestep forward for water simulation. 
          R,r      Reset positions and velocities for all particles.
        )";

    //about position and velocity
    //starting corner for the fluid particles;
    Particles particles;
    Eigen::Vector3d corner(-1., -1., -1.);
    Eigen::Vector3i num_points(20, 30, 20);
    double step_size = 0.054;
    const Eigen::RowVector3d particle_color(0.1, 0.9, 0.9);

    init_position(particles, corner, num_points, step_size);
    const auto& reset = [&]() {
        init_position(particles, corner, num_points, step_size);
        viewer.data_list[xid].set_points(particles.position, (1. - (1. - particle_color.array()) * .9));
    };

    const auto& move = [&]() {
        // TODO: fill the algorithm
        simulation_step(particles, double(0.008));
        viewer.data_list[xid].set_points(particles.position, (1. - (1. - particle_color.array()) * .9));
    };

    viewer.callback_key_pressed =
        [&](igl::opengl::glfw::Viewer&, unsigned char key, int) -> bool {
        switch (key) {
            case 'A':
            case 'a':
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

    viewer.data_list[xid].set_colors(particle_color);
    viewer.data_list[xid].set_points(particles.position, (1. - (1. - particle_color.array()) * .9));
    viewer.data_list[xid].point_size = 6.0;
    viewer.launch();

    

    return EXIT_SUCCESS;
}


