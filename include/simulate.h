#pragma once
#include "sim_data.h"
#include <Eigen/Sparse>
void simulate(sim_data& data, double dt);
void print_data(sim_data& data);
struct Vector3iCompare{
	bool operator() (const Eigen::Vector3i& lhs, const Eigen::Vector3i& rhs) const
	{
		return std::make_tuple(lhs(0), lhs(1), lhs(2)) < std::make_tuple(rhs(0), rhs(1), rhs(2));
	}
};

void set_grid_with_particles(sim_data& d);
void set_particles_with_grid(sim_data& d);

// calculate d vector: gridW*gridH*gridL (divergence of velocity at each point)
void calc_velocity_div(sim_data& data, Eigen::VectorXd& div);


// D * p should give vector of pressure gradients
// Pressure vector: p[x + y * gridW + z * gridW * gridH] = pressure at grid point x,y,z 
void calc_D_matrix(int gw, int gh, int gl, Eigen::SparseMatrix<double>& D);


// Build matrix B for LHS
// B * D * p = divergence of pressure gradient at each point
void calc_B_matrix(int gw, int gh, int gl, Eigen::SparseMatrix<double>& B);
