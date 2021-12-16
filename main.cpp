#include <init_particles.h>
#include <SPH_step.h>
#include <collision_response.h>

#include <igl/list_to_matrix.h>
#include <igl/opengl/glfw/Viewer.h>
#include <Eigen/Core>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cstdlib>

void fix_points(Eigen::MatrixXd& Position,
	Eigen::MatrixXd& Velocity,
	const Eigen::RowVector3d& grid_corner,
	const int nx,
	const int ny,
	const int nz,
	const double h)
{
	for (int i = 0; i < Position.rows(); i++) {
		if (Position(i, 0) < grid_corner(0)) {
			Position(i, 0) = grid_corner(0) + 1e-10;;
		}
		if (Position(i, 0) > grid_corner(0) + nx * h) {
			Position(i, 0) = grid_corner(0) + nx * h - 1e-10;

		}
		if (Position(i, 1) < grid_corner(1)) {
			Position(i, 1) = grid_corner(1) + 1e-10;

		}
		if (Position(i, 1) > grid_corner(1) + ny * h) {
			Position(i, 1) = grid_corner(1) + ny * h - 1e-10;

		}
		if (Position(i, 2) < grid_corner(2)) {
			Position(i, 2) = grid_corner(2) + 1e-10;

		}
		if (Position(i, 2) > grid_corner(2) + nz * h) {
			Position(i, 2) = grid_corner(2) + nz * h - 1e-10;
		}
	}

}
int main(int argc, char* argv[])
{

	Eigen::Vector3i num_points(15, 25, 15);
	double step_size = 0.05;
	Eigen::Vector3d corner(-0.5, -0.5, -0.5);
	const Eigen::RowVector3d particle_color(0.2, 0.8, 0.8);

	if (argc == 2 && argv[1] == std::string("SPH")) {
		igl::opengl::glfw::Viewer viewer;
		const int xid = viewer.selected_data_index;
		viewer.append_mesh();
        std::cout<<R"(
          A,a      Move one timestep forward for water simulation. 
          R,r      Reset positions and velocities for all particles.
        )";
		//about position and velocity
		//starting corner for the fluid particles;
		Particles particles;
		Coef coef;
		coef.MASS = 0.1;
		coef.STIFFNESS = 2.5;
		coef.RHO_IDEAL = 1000.;
		coef.VISCOSITY = 2.5;
		coef.H = 0.0635;
		coef.G = -9.8;
		coef.SUFRACE_TENSION = 1.0975;
		coef.SUFRACE_THRESH = 0.01;


		Eigen::MatrixXd N(6, 3);
		N <<
			0., -1., 0.,
			0., 1., 0.,
			1., 0., 0.,
			0., 0., 1.,
			-1., 0., 0.,
			0., 0., -1.;
		Eigen::MatrixXd P(6, 3);
		P <<
			0., 2., 0.,
			0., -1., 0.,
			-1., 0., 0.,
			0, 0., -1.,
			2., 0., 0.,
			0., 0., 2.;


		init_particles(particles, corner, num_points, step_size);
		const auto& reset = [&]()
		{
			init_particles(particles, corner, num_points, step_size);
			viewer.data_list[xid].set_points(particles.position, (1. - (1. - particle_color.array()) * .9));
		};

		const auto& move = [&]()
		{
			collision_response(particles, P, N);
			SPH_step(particles, 0.02, coef);

			viewer.data_list[xid].set_points(particles.position, (1. - (1. - particle_color.array()) * .9));

		};

		viewer.callback_key_pressed =
			[&](igl::opengl::glfw::Viewer&, unsigned char key, int)->bool
		{
			switch (key)
			{
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
	}


	return EXIT_SUCCESS;
}


