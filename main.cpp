#include <igl/list_to_matrix.h>
#include <igl/per_face_normals.h>
#include <igl/opengl/glfw/Viewer.h>
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

// viewer window
igl::opengl::glfw::Viewer viewer;

// the boundary of box
// TODO: change color
void draw_boudary(Particles& particles){
    // Find the bounding box
    Eigen::Vector3d m = particles.position.colwise().minCoeff() - 0.5 * Eigen::VectorXd::Ones(3*1);
    Eigen::Vector3d M = particles.position.colwise().maxCoeff() + 0.5 * Eigen::VectorXd::Ones(3*1);

    // Corners of the bounding box
    Eigen::MatrixXd V_box(8, 3);
    V_box << m(0), m(1), m(2),
        M(0), m(1), m(2),
        M(0), M(1), m(2),
        m(0), M(1), m(2),
        m(0), m(1), M(2),
        M(0), m(1), M(2),
        M(0), M(1), M(2),
        m(0), M(1), M(2);

    // Edges of the bounding box
    viewer.data().add_points(V_box,Eigen::RowVector3d(1,0,0));
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

    Eigen::Vector3d m = particles.position.colwise().minCoeff() - 0.5 * Eigen::VectorXd::Ones(3*1);
    Eigen::Vector3d M = particles.position.colwise().maxCoeff() + 0.5 * Eigen::VectorXd::Ones(3*1);

    // Corners of the bounding box
    Eigen::MatrixXd V_wall;
    V_wall.resize(8, 3);
    V_wall << m(0), m(1), m(2),
        M(0), m(1), m(2),
        M(0), M(1), m(2),
        m(0), M(1), m(2),
        m(0), m(1), M(2),
        M(0), m(1), M(2),
        M(0), M(1), M(2),
        m(0), M(1), M(2);

    Eigen::MatrixXi F_wall(6, 3);
    F_wall << 0, 1, 2,
            4, 7, 6, 
            0, 3, 7,
            6, 2, 1,
            4, 5, 1,
            3, 2, 6;

    Eigen::MatrixXd N_wall;
    igl::per_face_normals(V_wall, F_wall, N_wall); 


    const auto& reset = [&]() {
        init_position(particles, corner, num_points, step_size);
        viewer.data_list[xid].set_points(particles.position, (1. - (1. - particle_color.array()) * .9));
    };

    const auto& move = [&]() {
        // TODO: fill the algorithm
        simulation_step( particles, V_wall, F_wall, N_wall, double(0.008));
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


    draw_boudary(particles);
    viewer.data_list[xid].set_colors(particle_color);
    viewer.data_list[xid].set_points(particles.position, (1. - (1. - particle_color.array()) * .9));
    viewer.data_list[xid].point_size = 6.0;
    viewer.launch();

    return EXIT_SUCCESS;
}
