#include <collision_plane.h>

void collision_plane(
    Particles& particles,
    Eigen::MatrixXd& pred_position,
    const Eigen::MatrixXd &V_wall,
    const Eigen::MatrixXi &F_wall,
	const Eigen::MatrixXd &N_wall)
{

    int N = pred_position.rows();
    int num_walls = N_wall.rows();

    for(int i=0; i < N; i++){
        Eigen::Vector3d pos = pred_position.row(i);

        Eigen::Vector3d n;
        for(int w = 0; w < num_walls; w++){
            // n.dot(p - v)
            n = N_wall.row(w);
            Eigen::Vector3d v = V_wall.row(F_wall(w, 0));
            double diff = n.dot(pos - v);

            if(diff < 0){
                pos -= diff * n;
                pred_position.row(i) = pos;
            }

        }


    }

        

}