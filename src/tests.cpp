#include "tests.h"
#include <simulate.h>
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>
using namespace Eigen;
using namespace std;

sim_data initSimData(){
	sim_data data(3,3,3);
	cout << "initSimData\n";
	data.particles.push_back(particle(Vector3d(1.0,1.5,0.5), Vector3d(0,1,0)));
	//data.particles.push_back(particle(Vector3d(1.5,0.5,0.5), Vector3d(1.2,0,0)));
	
	//data.particles.push_back(particle(Vector3d(1.5,1.5,0.5), Vector3d(1,0,0)));
	//data.particles.push_back(particle(Vector3d(1.5,0.5,0.5), Vector3d(1.2,0,0)));
	//data.particles.push_back(particle(Vector3(2.5,1.5,0.5), Vector3d(1,0,0)));
	
	/*Vector3d origin = data.gridDims.cast<double>() * 0.5f;
	for(int i = 0; i < num; i++){
		//Vector3d offset = Vector3d::Zero();
		Vector3d offset = Vector3d::Random();
		offset(0) *= data.gridDims(0) * 0.25;
		offset(1) *= data.gridDims(1) * 0.25;
		offset(2) *= data.gridDims(2) * 0.25;
		data.particles.push_back(particle(origin + offset, Vector3d(0,1,0)));
		//data.particles.push_back(particle(origin + offset, Vector3d(0,1,0) + Vector3d::Random()));
	}*/
	return data;
}

void print_vector(VectorXd& v, int w, int h, int l){
	for(int z = 0; z < l; z ++){
		cout << "z = " << z << endl;
		for(int x = 0; x < w; x ++){
			for(int y = 0; y < h; y ++){
			
				cout << v(x + y * w + z * w * h) << " ";
			}
			cout << endl;
		}
		cout << endl;
	}
}

int idx(int x, int y, int z, Vector3i dims){
	return x + y * dims(0) + z * dims(0) * dims(1);
}

void runTests(){
	cout << "TEST START\n";	
	sim_data data = initSimData();
	int gw = data.gridDims(0);
	int gh = data.gridDims(1);
	int gl = data.gridDims(2);
	
	// Check vel distrib. on grid
	set_grid_with_particles(data);
	for(int z = 0; z <= data.gridDims(2); z++){
		cout << "gridVel[1][mid]:\n" << data.gridVel[1][z] << endl;
	}
	
	// Check redist.
	set_particles_with_grid(data);
	for(int i = 0; i < data.particles.size(); i++){
		cout << "part.vel: " << data.particles[i].vel.transpose() << endl;
	}
	
	// calculate d vector: gridW*gridH*gridL (divergence of velocity at each point)
	VectorXd div;
	calc_velocity_div(data, div);
	cout << "Divergence of velocity: \n";
	print_vector(div, gw, gh, gl);
	
	// D * p should give vector of pressure gradients
	// Pressure vector: p[x + y * gridW + z * gridW * gridH] = pressure at grid point x,y,z 
	VectorXd p = VectorXd::Zero(gw * gh * gl);
	Vector3i dims(gw, gh, gl);
	p(idx(1,1,1,dims)) = 1;
	
	SparseMatrix<double> D;
	calc_D_matrix(gw, gh, gl, D);
	//cout << "D matrix: \n" << MatrixXd(D);
	
	VectorXd pGrad = D * p;
	// print pressure gradients
	cout << "Pressure gradient: \n";
	for(int z = 0; z < gl; z ++){
		cout << "z = " << z << endl;
		for(int x = 0; x < gw; x ++){
			for(int y = 0; y < gh; y ++){
				int index = idx(x,y,z, dims) * 3;
				cout << "(" << pGrad(index) << "," <<  pGrad(index+1) << "," <<  pGrad(index+2) << ") ";
			}
			cout << endl;
		}
		cout << endl;
	}
	
	// Build matrix B for LHS
	// B * D * p = divergence of pressure gradient at each point
	SparseMatrix<double> B;
	calc_B_matrix(gw, gh, gl, B);
	VectorXd pDivergence = B * D * p;
	print_vector(pDivergence, gw, gh, gl);
	
	// Solve
	MatrixXd temp = B*D;
	SparseMatrix<double> A = temp.sparseView();
	
	LeastSquaresConjugateGradient<SparseMatrix<double>> solver;
	//LeastSquaresConjugateGradient<SparseMatrixd, Lower|Upper> solver;
	solver.compute(A);
	double density = 0.000001;
	double dt = 0.017;
	
	VectorXd b = div * density / dt; //test
	VectorXd pressure = solver.solve(b);
	if(solver.info() != Success){
		cout << "Solver failed! \n";
	}
	else{
		//cout << "SUCCESS! \n";
		//data.print_debug_once = true;
		//
	}
	
	cout << "TEST END\n";
}











