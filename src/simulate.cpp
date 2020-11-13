#include "simulate.h"
#include <Eigen/Core>
#include <algorithm>
#include <math.h>
#include <iostream>
#include <Eigen/Dense>

#include <Eigen/IterativeLinearSolvers>
//#include <EigenTypes.h>
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
	double w = (1-diff(0)) * (1-diff(1)) * (1-diff(2)); //?
	if (w < 0 || w > 1){
		//cout << "OUT OF RANGE WEIGHT!!!"<< endl;
	}
	return w;
}

void step_advection(sim_data& d, double dt){
	bool alert = false;
	for(int i = 0; i < d.particles.size(); i++){
		d.particles[i].pos += d.particles[i].vel * dt;
		Vector3d pos = d.particles[i].pos;
		if(pos(0) < 0 || pos(1) < 0 || pos(2) < 0){
			alert = true;
		}
	}
	if(alert){
		cout << "Particle < 0 !!\n";
	}
}

void step_external_force(sim_data& d, double dt){
	for(int i = 0; i < d.particles.size(); i++){
		d.particles[i].vel += d.gravity * dt;
	}
}


#include <map>
void set_grid_with_particles(sim_data& d){
	// May need to change!
	// check piazza "Bilinear Interpolation"
	// Consider: if 10 particles are near a point and moving in same direction
	// , does this mean that the particles would move essentially 10 times faster??
	//   => isn't normalization needed? (mentioned in comment, suggests to divide by sum of weights?)
	
	// use map
	map<Vector3i, double, Vector3iCompare> wSums;
	
	// zero out velocity grid!
	zero_out_grid(d.gridVel[0]);
	zero_out_grid(d.gridVel[1]);
	zero_out_grid(d.gridVel[2]);

	for(int i = 0; i < d.particles.size(); i++){typedef Triplet<double>td;
		particle part = d.particles[i];
		Vector3i iPos = clampI(part.pos, d.gridDims - Vector3i(1,1,1));
		Vector3d dPos = clampD(part.pos, d.gridDims.cast<double>());
		//cout << "iPos " << iPos << endl;
		//cout << "dPos " << dPos << endl;
		// loop over 8 points on grid to contribute to 
		for(int ox = 0; ox <= 1; ox ++){
			for(int oy = 0; oy <=1; oy ++){
				for(int oz = 0; oz <= 1; oz ++){
					Vector3i gPos = iPos + Vector3i(ox, oy, oz);
					double w = calc_weight(dPos, gPos.cast<double>());
					if(w != 0){
						if(wSums.count(gPos) == 0){
							wSums.insert(pair<Vector3i, double>(gPos, 0.0));
						}
						wSums[gPos] += w;
						//wSums[gPos] += w;
						for(int dim = 0; dim < 3; dim++){
							d.gridVel[dim][gPos[2]](gPos(0), gPos(1)) += w * part.vel[dim];
						}
					}
				}
			}
		}
	}
	
	if(d.normalize_velocity_grid){
		for(auto i = wSums.begin(); i != wSums.end(); ++i){
			for(int dim = 0; dim < 3; dim++){
				double weight = (*i).second;
				if(weight != 0){
					Vector3i vi = (*i).first;
					d.gridVel[dim][vi(2)](vi(0), vi(1)) /= (*i).second;
				}
			}
		}
	}
	
	if(d.print_debug_once){
		for(int z = 0; z <= d.gridDims(2); z++){
			cout << "gridVel[1][mid]:" << d.gridVel[1][z] << endl;
		}
		for(auto i = wSums.begin(); i != wSums.end(); ++i){
			cout << "wSum: " << (*i).second << endl;
		}
	}
}



void set_particles_with_grid(sim_data& d){
	for(int i = 0; i < d.particles.size(); i++){
		particle part = d.particles[i];
		d.particles[i].vel *= 0; // Zero out velocity
		Vector3i iPos = clampI(part.pos, d.gridDims - Vector3i(1,1,1));
		Vector3d dPos = clampD(part.pos, d.gridDims.cast<double>());
		if(d.print_debug_once){
			cout << "\npart.pos = " << part.pos.transpose() << endl;
			cout << "iPos = " << iPos.transpose() << endl;
			cout << "dPos = " << dPos.transpose() << endl;
			cout << "weights = ";
		}
		// loop over 8 points on grid to contribute to 
		for(int ox = 0; ox <= 1; ox ++){
			for(int oy = 0; oy <=1; oy ++){
				for(int oz = 0; oz <= 1; oz ++){
					Vector3i gPos = iPos + Vector3i(ox, oy, oz);
					double w = calc_weight(dPos, gPos.cast<double>());
					if(d.print_debug_once){
						cout << w << ", ";
					}
					for(int dim = 0; dim < 3; dim++){
						//d.gridVel[dim][gPos[2]](gPos(0), gPos(1)) += w * part.vel[dim];
						// DON'T use part, is copy of value
						//d.particles[i].vel[dim] += w * d.gridVel[dim][gPos[2]](gPos(0), gPos(1));
						d.particles[i].vel(dim) += w * d.gridVel[dim][gPos[2]](gPos(0), gPos(1));
					}
				}
			}
		}
		if(d.print_debug_once){
			cout << "Vel: " << d.particles[i].vel.transpose() << endl;
		}
	}
}

int getIndex(int x, int y, int z, Vector3i dims){
	return x + y * dims(0) + z * dims(0) * dims(1);
}

float getVel(int i, int x, int y, int z, std::vector<std::vector<Eigen::MatrixXd>>& gridVel){
	return gridVel[i][z](x,y);
}

float getPressure(int x, int y, int z, std::vector<Eigen::MatrixXd>& gridP){
	return gridP[z](x,y);
}

// calculate d vector: gridW*gridH*gridL (divergence of velocity at each point)
void calc_velocity_div(sim_data& data, VectorXd& div){
	int gw = data.gridDims(0);
	int gh = data.gridDims(1);
	int gl = data.gridDims(2);
	div = VectorXd::Zero(gw * gh * gl);
	int partCount = 0;
	for(int z = 0; z < gl; z++){
		for(int y = 0; y < gh; y++){
			for(int x = 0; x < gw; x++){
				// zero out some velocities in calculations if near edges for boundary condition
				/*double vx1 =  getVel(0,x,y,z,data.gridVel);
				double vx2 =  getVel(0,x+1,y,z,data.gridVel);
				double vy1 = getVel(1,x,y,z,data.gridVel);
				double vy2 = getVel(1,x,y+1,z,data.gridVel);
				double vz1 = getVel(2,x,y,z,data.gridVel);
				double vz2 = getVel(2,x,y,z+1,data.gridVel);*/
				
				double vx1 = x == 0 ? 0 : getVel(0,x,y,z,data.gridVel);
				double vx2 = x == gw-1 ? 0 : getVel(0,x+1,y,z,data.gridVel);
				double vy1 = y == 0 ? 0 : getVel(1,x,y,z,data.gridVel);
				double vy2 = y == gh-1 ? 0 : getVel(1,x,y+1,z,data.gridVel);
				double vz1 = z == 0 ? 0 : getVel(2,x,y,z,data.gridVel);
				double vz2 = z == gl-1 ? 0 : getVel(2,x,y,z+1,data.gridVel);
				
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
				/*int ax = tx == 0? 1 : tx; // adjusted for x gradient at edge
				Dlist.push_back(T(DRow, getIndex(ax-1,ty,tz,dims), -1));
				Dlist.push_back(T(DRow, getIndex(ax,ty,tz,dims), 1));
				DRow ++;
				
				int ay = ty == 0? 1 : ty;
				Dlist.push_back(T(DRow, getIndex(tx,ay-1,tz,dims), -1));
				Dlist.push_back(T(DRow, getIndex(tx,ay,tz,dims), 1));
				DRow ++;
				
				int az = tz == 0? 1 : tz;
				Dlist.push_back(T(DRow, getIndex(tx,ty,az-1,dims), -1));
				Dlist.push_back(T(DRow, getIndex(tx,ty,az,dims), 1));
				DRow ++;*/
				/*
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
				DRow ++;*/
			}
		}
	}
	D.setFromTriplets(Dlist.begin(), Dlist.end());
}

// Build matrix B for LHS
// B * D * p = divergence of pressure gradient at each point
void calc_B_matrix(int gw, int gh, int gl, SparseMatrixd& B){
	vector<T> Blist;
	B.resize(gw*gh*gl, 3*gw*gh*gl);
	Vector3i dims(gw, gh, gl);
	int BRow = 0;
	for(int tz = 0; tz < gl; tz++){
		for(int ty = 0; ty < gh; ty++){
			for(int tx = 0; tx < gw; tx++){
				// pick the first component of the current position's pressure gradient
				int ax = tx == 0 ? 1 : tx;
				if(tx != 0){
					Blist.push_back(T(BRow, getIndex(ax-1, ty, tz, dims)*3+0, -1)); 
				}
				if(tx != gw -1){
					Blist.push_back(T(BRow, getIndex(ax,ty,tz,dims)*3+0, 1));
				}
				
				int ay = ty == 0 ? 1 : ty;
				if(ty != 0){
					Blist.push_back(T(BRow, getIndex(tx, ay-1, tz, dims)*3+1, -1));
				}
				if(ty != gh - 1){
					Blist.push_back(T(BRow, getIndex(tx, ay, tz, dims)*3+1, 1));
				}
				
				int az = tz == 0 ? 1 : tz;
				if(tz != 0){
					Blist.push_back(T(BRow, getIndex(tx, ty, az-1, dims)*3+2, -1));
				}
				if(tz != gl - 1){
					Blist.push_back(T(BRow, getIndex(tx, ty, az, dims)*3+2, 1));
				}
				BRow ++;
				/*int ax = tx == gw - 1 ? gw -2 : tx;
				Blist.push_back(T(BRow, getIndex(ax, ty, tz, dims)*3+0, -1)); 
				Blist.push_back(T(BRow, getIndex(ax+1,ty,tz,dims)*3+0, 1));
				
				int ay = ty == gh - 1 ? gh- 2: ty;
				Blist.push_back(T(BRow, getIndex(tx, ay, tz, dims)*3+1, -1));
				Blist.push_back(T(BRow, getIndex(tx, ay+1, tz, dims)*3+1, 1));
				
				int az = tz == gl - 1 ? gl - 2 : tz;
				Blist.push_back(T(BRow, getIndex(tx, ty, az, dims)*3+2, -1));
				Blist.push_back(T(BRow, getIndex(tx, ty, az+1, dims)*3+2, 1));
				BRow ++;*/
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
	VectorXd div = VectorXd::Zero(gw * gh * gl);
	calc_velocity_div(data, div);
	
	// D * p should give vector of pressure gradients
	// grad[3*(x + y * gridW + z * gridW * gridH)] = gradient at x,y,z
	// Build 3*gridW*gridH*gridL X gridW*gridH*gridL matrix D for calcing gradient
	SparseMatrixd D(3*gw*gh*gl, gw*gh*gl);
	calc_D_matrix(gw, gh, gl, D);
	
	
	// Build matrix B for LHS
	// B * D * p = divergence of pressure gradient at each point
	SparseMatrixd B(gw*gh*gl, 3*gw*gh*gl);
	calc_B_matrix(gw,gh,gl, B);
	
	// B*D
	MatrixXd temp = B*D;
	SparseMatrixd A = temp.sparseView();
	
	LeastSquaresConjugateGradient<SparseMatrixd> solver;
	//ConjugateGradient<SparseMatrixd, Lower|Upper> solver;
	solver.compute(A);
	double density = 0.01;
	VectorXd b = div * density / dt;
	
	//cout << "b: " << b << endl;
	
	VectorXd pressure = solver.solve(b);
	
	if(temp.hasNaN()){
		cout << "A has NaN!!\n";
	}
	if(b.hasNaN()){
		cout << "b has NaN!!\n";
	}
	
	if(solver.info() != Success){
		cout << "Solver failed! \n";
		
	}
	else{
		//cout << "SUCCESS! \n";
		//data.print_debug_once = true;
		//
	}
	if(data.print_debug_once){
		cout << "Pressure: " << pressure.transpose() << endl;
	}
	
	// Need to put velocity update in
	// v = v - dt / density * grad(pressure)
	for(int tz = 0; tz < data.gridDims(2); tz++){
		for(int ty = 0; ty < data.gridDims(1); ty++){
			for(int tx = 0; tx < data.gridDims(0); tx++){
			// TODO: fix edge issue
			    int ax = tx == gw -1 ? gw - 2 : tx; // adjusted for x gradient at edge
			    int ay = ty == gh -1 ? gh - 2 : ty; // adjusted for y gradient at edge
			    int az = tz == gl -1 ? gl - 2 : tz; // adjusted for z gradient at edge
				
				/*int x = min(data.gridDims(0)-2, tx); // gradient at the right/top/far edges is a repeated value
				int y = min(data.gridDims(1)-2, ty);
				int z = min(data.gridDims(2)-2, tz);*/
				
				// get pressure gradient
				Vector3d pGrad = Vector3d::Zero();
				
				pGrad(0) = pressure(getIndex(ax+1,ty,tz,data.gridDims)) - pressure(getIndex(ax,ty,tz,data.gridDims));
				pGrad(1) = pressure(getIndex(tx,ay+1,tz,data.gridDims)) - pressure(getIndex(tx,ay,tz,data.gridDims));
				pGrad(2) = pressure(getIndex(tx,ty,az+1,data.gridDims)) - pressure(getIndex(tx,ty,az,data.gridDims));
				
				
				data.gridVel[0][tz](tx,ty) -= dt / density * pGrad(0);
				data.gridVel[1][tz](tx,ty) -= dt / density * pGrad(1);
				data.gridVel[2][tz](tx,ty) -= dt / density * pGrad(2);
				
				
				if(false){//data.print_debug_once){
					cout << "-dt / density * pGrad: " << -dt / density * pGrad.transpose() << endl;
					Vector3d tempVel(data.gridVel[0][tz](tx,ty), data.gridVel[1][tz](tx,ty), data.gridVel[2][tz](tx,ty));
					cout << "data.gridVel: " << tempVel.transpose() << endl;;
				}
			}
		}
	}
}

void simulate(sim_data& data, double dt){
	step_advection(data, dt);
	step_external_force(data, dt);
	
	// set grid with velocities from particles
	set_grid_with_particles(data);
	
	// zero out pressure grid first?
	// solve for pressure, update gridVel
	pressure_projection(data, dt);
	
	set_particles_with_grid(data);
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


/*
	// Build q vector
	VectorXd q(3 * (data.gridDims(0) + 1)* (data.gridDims(1) + 1)* (data.gridDims(2) + 1));
	int partCount = 0;
	for(int z = 0; z <= data.gridDims(2); z++){
		for(int y = 0; y <= data.gridDims(1); y++){
			for(int x = 0; x <= data.gridDims(0); x++){
				int ind = getIndex(x,y,z, data.gridDims + Vector3i(1,1,1));
				q(ind*3) = data.particles[partCount].vel(0);
				q(ind*3+1) = data.particles[partCount].vel(1);
				q(ind*3+2) = data.particles[partCount].vel(2);
			}
		}
	}
	
	// Build RHS B Matrix
	vector<T> Blist;
	SparseMatrixd B(, q.rows());
	int count = 0;
	for(int z = 0; z < data.gridDims(2); z++){
		for(int y = 0; y < data.gridDims(1); y++){
			for(int x = 0; x < data.gridDims(0); x++){
				int curInd = 3 * getIndex(x, y, z, data.gridDims);
				Blist.push_back(T(count, curInd + 0, -1)); // x velocity  at current position
				Blist.push_back(T(count, 3 * getIndex(x+1, y, z, data.gridDims), 1));
				
				Blist.push_back(T(count, curInd + 1, -1)); // y velocity at current position
				Blist.push_back(T(count, 3 * getIndex(x,y+1,z,data.gridDims) + 1, 1));
				
				Blist.push_back(T(count, curInd + 2, -1)); // z velocity at current pos
				Blist.push_back(T(count, 3 * getIndex(x,y,z+1,data.gridDims) + 2, 1));
				count ++;
			}
		}
	}*/
