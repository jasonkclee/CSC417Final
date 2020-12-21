#pragma once
#include <vector>
#include <Eigen/Core>
#include "particle.h"

void build_grid(std::vector<Eigen::MatrixXd>& grid, int w, int h, int l);


struct sim_data{
	bool print_debug_once;
	bool flip_method;
	Eigen::Vector3d gravity;
	Eigen::Vector3i gridDims;
	std::vector<particle> particles;
	
	
	std::vector<std::vector<Eigen::MatrixXd>> gridVel, dGridVel;

	sim_data(int gridW, int gridH, int gridL): gravity(0,-0.1,0), print_debug_once(false),
		gridDims(gridW, gridH, gridL), flip_method(true){
		for(int d = 0; d < 3; d++){
			gridVel.push_back(std::vector<Eigen::MatrixXd>());
			dGridVel.push_back(std::vector<Eigen::MatrixXd>());
		}
		build_grid(gridVel[0], gridW + 1, gridH + 2, gridL + 2);
		build_grid(gridVel[1], gridW + 2, gridH + 1, gridL + 2);
		build_grid(gridVel[2], gridW + 2, gridH + 2, gridL + 1);
		
		build_grid(dGridVel[0], gridW + 1, gridH + 2, gridL + 2);
		build_grid(dGridVel[1], gridW + 2, gridH + 1, gridL + 2);
		build_grid(dGridVel[2], gridW + 2, gridH + 2, gridL + 1);
		
	}
	
	
};

typedef struct sim_data sim_data;
