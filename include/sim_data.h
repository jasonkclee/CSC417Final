#pragma once
#include <vector>
#include <Eigen/Core>
#include "particle.h"
struct sim_data{
	bool print_debug_once;
	bool normalize_velocity_grid;

	sim_data(int gridW, int gridH, int gridL, double gridScale): gravity(0,-1.0,0), gridScale(gridScale), print_debug_once(false),
		gridDims(gridW, gridH, gridL), normalize_velocity_grid(false){
		for(int d = 0; d < 3; d++){
			gridVel.push_back(std::vector<Eigen::MatrixXd>());
		}
		for(int i = 0; i < gridL; i++){
			gridPressure.push_back(Eigen::MatrixXd::Zero(gridW, gridH));
			for(int d = 0; d < 3; d++){
				gridVel[d].push_back(Eigen::MatrixXd::Zero(gridW+1, gridH+1));
			}
		}
		for(int d = 0; d < 3; d++){
				gridVel[d].push_back(Eigen::MatrixXd::Zero(gridW+1, gridH+1));
		}
	}
	Eigen::Vector3i gridDims;
	std::vector<particle> particles;
	Eigen::Vector3d gravity;
	
	double gridScale;
	std::vector<std::vector<Eigen::MatrixXd>> gridVel;
	std::vector<Eigen::MatrixXd> gridPressure;
	
};

typedef struct sim_data sim_data;
