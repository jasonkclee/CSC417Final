#pragma once
#include <Eigen/Core>
struct particle{
	particle() : pos(0,0,0), vel(0,0,0){}
	particle(Eigen::Vector3d p, Eigen::Vector3d v) : pos(p), vel(v){}
	Eigen::Vector3d pos, vel;
};
typedef struct particle particle;
