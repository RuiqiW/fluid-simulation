#include <collision_mesh.h>
#include <igl/signed_distance.h>

void collision_mesh(
    Particles& particles,
    Eigen::MatrixXd& pred_position,
    const Eigen::MatrixXd &V_obj,
    const Eigen::MatrixXi &F_obj,
    igl::AABB<Eigen::MatrixXd, 3> &tree
){
    Eigen::VectorXd sqrD;
    Eigen::VectorXi I;
    Eigen::MatrixXd C;

    tree.squared_distance(V_obj,F_obj, pred_position, sqrD,I,C);


    Eigen::VectorXd S_dist; // signed distance
    Eigen::MatrixXd N;

    igl::SignedDistanceType sign_type = igl::SIGNED_DISTANCE_TYPE_PSEUDONORMAL;
    igl::signed_distance(pred_position, V_obj, F_obj, sign_type, S_dist, I, C, N);


    for(int idx=0; idx < pred_position.rows(); idx++){
        if(S_dist(idx) < 0){
            Eigen::Vector3d pos = pred_position.row(idx);
            Eigen::Vector3d n = N.row(idx);
            pos -= S_dist(idx) * n;
            pred_position.row(idx) = pos;
        }
    }
}