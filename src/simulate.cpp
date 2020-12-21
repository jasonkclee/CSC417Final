#include "simulate.h"
#include <Eigen/Core>
#include <algorithm>
#include <math.h>
#include <iostream>
#include <Eigen/Dense>

#include <Eigen/IterativeLinearSolvers>
typedef Eigen::SparseMatrix<double>SparseMatrixd;
using namespace Eigen;
using namespace std;
typedef Triplet<double>T;

Vector3i clampI(Vector3d pos, Vector3i pMax){
    Vector3i pMin = Vector3i::Zero();
	Vector3i iPos((int)std::floor(pos(0)), (int)std::floor(pos(1)), (int)std::floor(pos(2)));
	iPos(0) = min(pMax(0), max(pMin(0), iPos(0)));
	iPos(1) = min(pMax(1), max(pMin(1), iPos(1)));
	iPos(2) = min(pMax(2), max(pMin(2), iPos(2)));
	return iPos;
}

Vector3d clampD(Vector3d pos, Vector3d pMax){
	Vector3d vPos;
    Vector3d pMin = Vector3d::Zero();
	vPos(0) = min(pMax(0), max(pMin(0), pos(0)));
	vPos(1) = min(pMax(1), max(pMin(1), pos(1)));
	vPos(2) = min(pMax(2), max(pMin(2), pos(2)));
	return vPos;
}

void zero_out_grid(vector<MatrixXd>& grid){
	for(int i = 0; i < grid.size(); i++){
		grid[i] *= 0;
	}
}

double calc_weight(Vector3d v1, Vector3d v2){
	Vector3d diff = (v1 - v2).cwiseAbs();
	double w = (1-diff(0)) * (1-diff(1)) * (1-diff(2));
	return w;
}

void remove_particles(sim_data& d){
	auto it = d.particles.begin();
	while(it != d.particles.end()){
		particle p = *it;
		if(p.pos(1) < 0.0){
			it = d.particles.erase(it);
		}
		else{
			++it;
		}
	}
}

void step_advection(sim_data& d, double dt){
	for(int i = 0; i < d.particles.size(); i++){
		d.particles[i].pos += d.particles[i].vel * dt;
		Vector3d pos = d.particles[i].pos;
	}
}

void step_external_force(sim_data& d, double dt){
	for(int i = 0; i < d.particles.size(); i++){
		d.particles[i].vel += d.gravity * dt;
	}
}

void set_grid_with_particles(sim_data& d){
	
	// zero out velocity grid!
	zero_out_grid(d.gridVel[0]);
	zero_out_grid(d.gridVel[1]);
	zero_out_grid(d.gridVel[2]);
	
	// need to consider different offset for every grid
	for(int dim = 0; dim < 3; dim++){
		Vector3d offset = Vector3d::Zero();
		int dim2 = (dim+1)%3;
		int dim3 = (dim+2)%3;
		offset(dim2) = 0.5;
		offset(dim3) = 0.5;
		for(int i = 0; i < d.particles.size(); i++){
			particle part = d.particles[i];
			// position offsetted
			Vector3d dposForIndex = part.pos + offset;
			Vector3i maxIndex = Vector3i::Zero();
			maxIndex(dim) = d.gridDims(dim)-1;
			maxIndex(dim2) = d.gridDims(dim2);
			maxIndex(dim3) = d.gridDims(dim3);
			Vector3i iposForIndex = clampI(dposForIndex, maxIndex);
			for(int ox = 0; ox <= 1; ox ++){
				for(int oy = 0; oy <=1; oy ++){
					for(int oz = 0; oz <= 1; oz ++){
					
						Vector3i vectorIndex = iposForIndex + Vector3i(ox, oy, oz);
						double w = calc_weight(dposForIndex, vectorIndex.cast<double>());
						d.gridVel[dim][vectorIndex[2]](vectorIndex(0), vectorIndex(1)) += w * part.vel(dim);
					}
				}
			}
		}
	}
	
}



void set_particles_with_grid(sim_data& d){
	for(int i = 0; i < d.particles.size(); i++){
		d.particles[i].vel *= 0; // Zero out velocity
	}
	for(int dim = 0; dim < 3; dim++){
		Vector3d offset = Vector3d::Zero();
		int dim2 = (dim+1)%3;
		int dim3 = (dim+2)%3;
		offset(dim2) = 0.5;
		offset(dim3) = 0.5;
		for(int i = 0; i < d.particles.size(); i++){
			particle part = d.particles[i];
			
			// position offsetted
			Vector3d dposForIndex = part.pos + offset;
			Vector3i maxIndex = Vector3i::Zero();
			maxIndex(dim) = d.gridDims(dim)-1;
			maxIndex(dim2) = d.gridDims(dim2);
			maxIndex(dim3) = d.gridDims(dim3);
			Vector3i iposForIndex = clampI(dposForIndex, maxIndex);
			for(int ox = 0; ox <= 1; ox ++){
				for(int oy = 0; oy <=1; oy ++){
					for(int oz = 0; oz <= 1; oz ++){
						Vector3i vectorIndex = iposForIndex + Vector3i(ox, oy, oz);
						double w = calc_weight(dposForIndex, vectorIndex.cast<double>());
						d.particles[i].vel(dim) += w * d.gridVel[dim][vectorIndex[2]](vectorIndex(0), vectorIndex(1));
					}
				}
			}
		}
	}
}

void set_particles_with_grid_flip(sim_data& d){
	
	for(int dim = 0; dim < 3; dim++){
		Vector3d offset = Vector3d::Zero();
		int dim2 = (dim+1)%3;
		int dim3 = (dim+2)%3;
		offset(dim2) = 0.5;
		offset(dim3) = 0.5;
		for(int i = 0; i < d.particles.size(); i++){
			particle part = d.particles[i];
			
			// position offsetted
			Vector3d dposForIndex = part.pos + offset;
			Vector3i maxIndex = Vector3i::Zero();
			maxIndex(dim) = d.gridDims(dim)-1;
			maxIndex(dim2) = d.gridDims(dim2);
			maxIndex(dim3) = d.gridDims(dim3);
			Vector3i iposForIndex = clampI(dposForIndex, maxIndex);
			for(int ox = 0; ox <= 1; ox ++){
				for(int oy = 0; oy <=1; oy ++){
					for(int oz = 0; oz <= 1; oz ++){
						Vector3i vectorIndex = iposForIndex + Vector3i(ox, oy, oz);
						double w = calc_weight(dposForIndex, vectorIndex.cast<double>());
						d.particles[i].vel(dim) += w * d.dGridVel[dim][vectorIndex[2]](vectorIndex(0), vectorIndex(1));
					}
				}
			}
		}
	}
}

int getIndex(int x, int y, int z, Vector3i dims){
	return x + y * dims(0) + z * dims(0) * dims(1);
}

float getVel(int i, int x, int y, int z, std::vector<std::vector<Eigen::MatrixXd>>& gridVel){
	return gridVel[i][z](x,y);
}


// calculate d vector: gridW*gridH*gridL (divergence of velocity at each point)
void calc_velocity_div(sim_data& data, VectorXd& div){
	int gw = data.gridDims(0);
	int gh = data.gridDims(1);
	int gl = data.gridDims(2);
	div = VectorXd::Zero((gw+2) * (gh+2) * (gl+2));
	int partCount = 0;
	for(int z = 0; z < gl+2; z++){
		for(int y = 0; y < gh+2; y++){
			for(int x = 0; x < gw+2; x++){
				// zero out some velocities in calculations if near edges for boundary condition
				
				double vx1 = x == 0 ? 0 : getVel(0,x-1,y,z,data.gridVel);
				double vx2 = x == gw+1 ? 0 : getVel(0,x,y,z,data.gridVel);
				double vy1 = y == 0 ? 0 : getVel(1,x,y-1,z,data.gridVel);
				double vy2 = y == gh+1 ? 0 : getVel(1,x,y,z,data.gridVel);
				double vz1 = z == 0 ? 0 : getVel(2,x,y,z-1,data.gridVel);
				double vz2 = z == gl+1 ? 0 : getVel(2,x,y,z,data.gridVel);
				
				div(partCount) = vx2 - vx1
							   + vy2 - vy1
							   + vz2 - vz1;
				partCount ++;
			}
		}
	}
}

// D * p should give vector of pressure gradients
// Pressure vector: p[x + y * gridW + z * gridW * gridH] = pressure at grid point x,y,z 
void calc_D_matrix(int gw, int gh, int gl, SparseMatrixd& D){
	Vector3i dims(gw, gh, gl);
	vector<T> Dlist;
	D.resize(3*gw*gh*gl, gw*gh*gl);
	int DRow = 0;
	for(int tz = 0; tz < gl; tz++){
		for(int ty = 0; ty < gh; ty++){
			for(int tx = 0; tx < gw; tx++){
				// don't zero out gradients here; perform boundary step at calc_B_matrix
				int ax = tx == gw -1 ? gw - 2 : tx; // adjusted for x gradient at edge
				Dlist.push_back(T(DRow, getIndex(ax,ty,tz,dims), -1));
				Dlist.push_back(T(DRow, getIndex(ax+1,ty,tz,dims), 1));
				
				DRow ++;
				
				int ay = ty == gh-1? gh - 2 : ty;
				Dlist.push_back(T(DRow, getIndex(tx,ay,tz,dims), -1));
				Dlist.push_back(T(DRow, getIndex(tx,ay+1,tz,dims), 1));
				DRow ++;
				
				int az = tz == gl-1? gl - 2 : tz;
				Dlist.push_back(T(DRow, getIndex(tx,ty,az,dims), -1));
				Dlist.push_back(T(DRow, getIndex(tx,ty,az+1,dims), 1));
				DRow ++;
			}
		}
	}
	D.setFromTriplets(Dlist.begin(), Dlist.end());
}

// Build matrix B for LHS
// B * D * p = divergence of pressure gradient at each point
void calc_B_matrix(int gw, int gh, int gl, SparseMatrixd& B){
	vector<T> Blist;
	B.resize((gw+2) * (gh+2) * (gl+2), 3*(gw+2) * (gh+2) * (gl+2));
	Vector3i dims((gw+2), (gh+2), (gl+2));
	int BRow = 0;
	for(int tz = 0; tz < gl+2; tz++){
		for(int ty = 0; ty < gh+2; ty++){
			for(int tx = 0; tx < gw+2; tx++){
				// zero out/ignore some of the pressure gradients must be zeroed out due to wall
			
				if(tx != 0){
					Blist.push_back(T(BRow, getIndex(tx-1, ty, tz, dims)*3+0, -1));
				}
				if(tx != gw +1){
					Blist.push_back(T(BRow, getIndex(tx,ty,tz,dims)*3+0, 1));
				}
				
				if(ty != 0){
					Blist.push_back(T(BRow, getIndex(tx, ty-1, tz, dims)*3+1, -1));
				}
				if(ty != gh + 1){
					Blist.push_back(T(BRow, getIndex(tx, ty, tz, dims)*3+1, 1));
				}
				
				if(tz != 0){
					Blist.push_back(T(BRow, getIndex(tx, ty, tz-1, dims)*3+2, -1));
				}
				if(tz != gl + 1){
					Blist.push_back(T(BRow, getIndex(tx, ty, tz, dims)*3+2, 1));
				}
				BRow ++;
			}
		}
	}
	B.setFromTriplets(Blist.begin(), Blist.end());
	

}

void pressure_projection(sim_data& data, double dt){
	int gw = data.gridDims(0);
	int gh = data.gridDims(1);
	int gl = data.gridDims(2);
	
	// Pressure vector: pressure[x + y * gridW + z * gridW * gridH] = pressure at grid point x,y,z 
	// velocity vector q: velocity[3*(x + y * gridW + z * gridW * gridH) + i] = vel(i) at grid point x,y,z
	
	// calculate d vector: gridW*gridH*gridL (divergence of velocity at each point)
	VectorXd div;// = VectorXd::Zero(gw * gh * gl);
	calc_velocity_div(data, div);
	
	// D * p should give vector of pressure gradients
	// grad[3*(x + y * gridW + z * gridW * gridH)] = gradient at x,y,z
	// Build 3*gridW*gridH*gridL X gridW*gridH*gridL matrix D for calcing gradient
	SparseMatrixd D;//(3*gw*gh*gl, gw*gh*gl);
	calc_D_matrix(gw, gh, gl, D);
	
	
	// Build matrix B for LHS
	// B * D * p = divergence of pressure gradient at each point
	SparseMatrixd B;//(gw*gh*gl, 3*gw*gh*gl);
	calc_B_matrix(gw,gh,gl, B);
	
	// B*D
	//MatrixXd temp = B*D;
	MatrixXd temp = B;
	SparseMatrixd A = temp.sparseView();
	
	LeastSquaresConjugateGradient<SparseMatrixd> solver;
	solver.compute(A);
	double density = 0.01;
	VectorXd b = div * density / dt;
	
	VectorXd pressure = solver.solve(b);
	
	if(solver.info() != Success){
		//cout << "Solver failed! \n";
	}
	
	if(data.print_debug_once){
		cout << "Pressure: " << pressure.transpose() << endl;
	}
	
	// Need to put velocity update in
	// v = v - dt / density * grad(pressure)
	zero_out_grid(data.dGridVel[0]);
	zero_out_grid(data.dGridVel[1]);
	zero_out_grid(data.dGridVel[2]);
	
	for(int dim = 0; dim < 3; dim++){
		int dim2 = (dim+1)%3;
		int dim3 = (dim+2)%3;
		for(int iz = 0; iz < data.gridVel[dim].size(); iz++){
			for(int iy = 0; iy < data.gridVel[dim][0].cols(); iy ++){
				for(int ix = 0; ix < data.gridVel[dim][0].rows(); ix ++){
					Vector3i idx(ix,iy,iz);
					Vector3i pressureDim((gw+2), (gh+2), (gl+2));
					Vector3i tIdx = idx;
					//tIdx(dim) += 1;
					
					//double pGrad = pressure(getIndex(ix,iy,iz,pressureDim)*3 + dim);
					double pGrad = pressure(getIndex(tIdx(0),tIdx(1),tIdx(2),pressureDim)*3 + dim);
					data.gridVel[dim][iz](ix,iy) -= dt / density * pGrad;
					data.dGridVel[dim][iz](ix,iy) -= dt / density * pGrad;
					
				}
			}
		}
	}
}

void simulate(sim_data& data, double dt){
	//remove_particles(data);
	step_advection(data, dt);
	step_external_force(data, dt);
	
	// set grid with velocities from particles
	set_grid_with_particles(data);
	
	// solve for pressure, update gridVel
	pressure_projection(data, dt);
	
	if(data.flip_method){
		set_particles_with_grid_flip(data);
	}
	else{
		set_particles_with_grid(data);
	}
	data.print_debug_once = false;
}

void print_data(sim_data& d){
	for(int z = 0; z <= d.gridDims(2); z++){
		cout << "gridVel[1][mid]:" << d.gridVel[1][z] << endl;
	}
	for(int i = 0; i < d.particles.size(); i++){
		particle part = d.particles[i];
		cout << "vel: " << part.vel.transpose() << endl;
	}
	
}


