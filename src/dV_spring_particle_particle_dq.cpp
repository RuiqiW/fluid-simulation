#include <dV_spring_particle_particle_dq.h>

void dV_spring_particle_particle_dq(Eigen::Ref<Eigen::Vector6d> f, Eigen::Ref<const Eigen::Vector3d> q0,  Eigen::Ref<const Eigen::Vector3d>     q1, double l0, double stiffness) {
    double len = (q1 - q0).norm();  
    Eigen::Vector6d q;  
    q << q0, q1;  
  
    Eigen::MatrixXd B(3, 6);  
    B << -1, 0, 0, 1, 0, 0,  
         0, -1, 0, 0, 1, 0,  
         0, 0, -1, 0, 0, 1;  
  
    // k * (qTBTBq^(1/2) - l0) * (qTBTBq)^(-1/2) * BTBq  
    f = stiffness * (len - l0) / len * B.transpose() * B * q;  
   
    
}