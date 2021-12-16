#ifndef ASSIGNMENT_SETUP_H
#define ASSIGNMENT_SETUP_H
//A6 code
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <EigenTypes.h>

#include <visualization.h>
#include <init_state_rigid_bodies.h>
#include <dV_spring_particle_particle_dq.h>

#define RB_OFFSET 12
#define RB_POS_OFFSET 9

//rigid body geometry
std::vector<std::pair<Eigen::MatrixXd, Eigen::MatrixXi>> geometry;

Eigen::VectorXd *q_ptr, *qdot_ptr;

//rigid body mass matrices
std::vector<Eigen::Matrix66d> mass_matrices;

//skinning matrix for rigid body (necesary for the way I interface with the libigl viewer)
Eigen::SparseMatrixd N;

double density = 1.0;

//selection spring
double k_selected = 1;

//some forces for the rigid bodies
Eigen::VectorXd forces; 

//a reminder that rigid body state is stored as 3x3 rotation matrix (body to world) and a 3x1 vector (world space posirion of center of mass)
//velocity is stored as a 3x1 angular velocity vector in world space and a 3x1 linear velocity vector in world space 
inline void simulate(Eigen::VectorXd &q, Eigen::VectorXd &qdot, double dt, double t) {  

    //Interaction spring
    Eigen::Vector3d mouse;
    Eigen::Vector6d dV_mouse;
    Eigen::Matrix36d rb_jacobian;

    forces.setZero();

    unsigned int irb = 0;
    for(unsigned int pickedi = 0; pickedi < Visualize::picked_vertices().size(); pickedi++) {   
        Eigen::Matrix3d R = Eigen::Map<const Eigen::Matrix3d>(q.segment<9>(12*irb).data());
        Eigen::Vector3d p = Eigen::Map<const Eigen::Vector3d>(q.segment<3>(12*irb + 9).data());

        mouse = Visualize::geometry(irb).row(Visualize::picked_vertices()[pickedi]).transpose() + Visualize::mouse_drag_world() + Eigen::Vector3d::Constant(1e-6);
        dV_spring_particle_particle_dq(dV_mouse, mouse, Visualize::geometry(irb).row(Visualize::picked_vertices()[pickedi]).transpose(), 0.0, (Visualize::is_mouse_dragging() ? k_selected : 0.));
        
        //std::cout<<dV_mouse.transpose()<<"\n";

        // rigid_body_jacobian(rb_jacobian, R, p, Visualize::geometry(irb).row(Visualize::picked_vertices()[pickedi]).transpose());
        // forces.segment<6>(6*irb) -= rb_jacobian.transpose()*dV_mouse.segment<3>(3);
    }



    // exponential_euler(q, qdot, dt, mass_matrices, forces);
}

inline void draw(Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::VectorXd> qdot, double t) {

    //update vertex positions using simulation
    for(unsigned int irb=0; irb<geometry.size();++irb) {
        Eigen::Map<const Eigen::Matrix3d> R = Eigen::Map<const Eigen::Matrix3d>(q.segment<9>(RB_OFFSET*irb).data());
        Eigen::MatrixXd V_rb = R*(geometry[irb].first.transpose());
        
        for(unsigned int jj=0; jj<V_rb.cols(); ++jj) {
            V_rb.col(jj) += q.segment<3>(RB_OFFSET*irb + RB_POS_OFFSET);
        }
        Visualize::update_vertex_positions(irb, Eigen::Map<const Eigen::VectorXd>(V_rb.data(), 3*V_rb.cols()));
    }

}

bool key_down_callback(igl::opengl::glfw::Viewer &viewer, unsigned char key, int modifiers) {

    if(key == 'R') {
        //reset the simulation
        forces.setZero();
        init_state_rigid_bodies(*q_ptr, *qdot_ptr,geometry.size());
        qdot_ptr->segment<3>(0) << 0.1, 0.1, 0.0;
    }
    return false;
}
inline void assignment_setup(int argc, char **argv, Eigen::VectorXd &q, Eigen::VectorXd &qdot) {

    q_ptr = &q;
    qdot_ptr = &qdot;

    //load geometric data 
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    igl::readOBJ("../data/torus.obj", V, F);

    geometry.push_back(std::make_pair(V,F));

    //some forces to show the integrator works 
    forces.resize(6*geometry.size());
    forces.setZero();

    //setup simulation 
    init_state_rigid_bodies(q,qdot,geometry.size());

    Eigen::Matrix3d inertia;
    Eigen::Matrix66d mass_matrix;
    double mass = 0;
    Eigen::Vector3d com;
    Eigen::Vector3d gravity;
    gravity << 0., -0.0000098, 0.;

    // for(unsigned int irb=0; irb < geometry.size(); ++irb) {
    //     inertia_matrix(inertia, com, mass, geometry[irb].first, geometry[irb].second, density);
        
    //     mass_matrix.setZero();
    //     mass_matrix.block(0,0,3,3) = inertia;
    //     mass_matrix(3,3) = mass;
    //     mass_matrix(4,4) = mass;
    //     mass_matrix(5,5) = mass;

    //     mass_matrices.push_back(mass_matrix);
        
    //     //inital angular velocity
    //     qdot.segment<3>(6*irb) << 0.1, 0.1, 0.0;

    //     //setup rigid bodies initial position in space 
    //     q.segment<3>(RB_OFFSET*irb + RB_POS_OFFSET) = com;

    //     forces.segment<3>(6*irb + 3) = mass_matrices[irb].block(3,3,3,3)*gravity;
    //     N.resize(geometry[irb].first.rows(), geometry[irb].first.rows());
    //     N.setIdentity();
    //     Visualize::add_object_to_scene(geometry[irb].first, geometry[irb].second, geometry[irb].first, geometry[irb].second, N,Eigen::RowVector3d(244,165,130)/255.);

    //     //fix up mesh
    //     for(unsigned int jj=0; jj<geometry[irb].first.rows(); ++jj) {
    //         geometry[irb].first.row(jj) -= com.transpose();
    //     }
    // }

    Visualize::viewer().callback_key_down = key_down_callback;

}

#endif

