#include <init_boundaries.h>
#include <igl/readOBJ.h>


void init_rabbit(
	Rabbit_Mesh &rabbit,
    const char* objFile,
    bool isUsed
) {
    // bool isUsed = true;
    if (isUsed == false){
        rabbit.isUsed = false;
        return;
    }
    // Eigen::MatrixXd V_obj;
    // Eigen::MatrixXi F_obj;
    // string objFile = "../data/coarser_bunny.obj";
    igl::readOBJ(objFile, rabbit.V_obj, rabbit.F_obj);
    rabbit.V_obj *= 0.005;
    Eigen::VectorXd tmp;
    tmp.resize(rabbit.V_obj.rows());
    tmp.setOnes();
    tmp *= 2.2;
    rabbit.V_obj.col(1) -= tmp;
    // rabbit.V_obj = V_obj;
    // rabbit.F_obj = F_obj;
    rabbit.tree.init(rabbit.V_obj,rabbit.F_obj);
    rabbit.isUsed = true;



}


void init_wall(
	Wall_Plane &wall,
    Eigen::Vector3d &m,
    Eigen::Vector3d &M,
    bool isUsed
) {
    if (isUsed == false){
        wall.isUsed = false;
        return;
    }
    wall.isUsed = true;

    
    wall.V_wall.resize(8, 3);
    wall.V_wall << m(0), m(1), m(2),
        M(0), m(1), m(2),
        M(0), M(1), m(2),
        m(0), M(1), m(2),
        m(0), m(1), M(2),
        M(0), m(1), M(2),
        M(0), M(1), M(2),
        m(0), M(1), M(2);

    wall.F_wall.resize(6,3);
    wall.F_wall << 0, 1, 2,
            4, 7, 6, 
            0, 3, 7,
            6, 2, 1,
            4, 5, 1,
            3, 2, 6;
}
