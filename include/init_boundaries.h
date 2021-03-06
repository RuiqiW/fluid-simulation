#include <boundary_property.h>
#include <Eigen/Dense>




void init_rabbit(
	Rabbit_Mesh &rabbit,
    const char* objFile,
    bool isUsed // whether the wall is in the scene
);


void init_wall(
	Wall_Plane &wall,
    Eigen::Vector3d & m, // min_boundary,
    Eigen::Vector3d &M, // max_boundary,
    bool isUsed // whether the wall is in the scene
);