#include <igl/opengl/glfw/Viewer.h>
#include <iostream>
#include "particle.h"
#include "sim_data.h"
#include "simulate.h"
#include "tests.h"
using namespace Eigen;
using namespace std;

MatrixXd getCubeV(double w, Vector3d center){
	MatrixXd V(8,3);
	V << 0.0,0.0,0.0,
    0.0,0.0,w,
    0.0,w,0.0,
    0.0,w,w,
    w,0.0,0.0,
    w,0.0,w,
    w,w,0.0,
    w,w,w;
    V.col(0) += VectorXd::Ones(V.rows()) * (center(0)-w * 0.5);
    V.col(1) += VectorXd::Ones(V.rows()) * (center(1)-w * 0.5);
    V.col(2) += VectorXd::Ones(V.rows()) * (center(2)-w * 0.5);
    return V;
}

MatrixXi getCubeF(){
    return (Eigen::MatrixXi(12,3)<<
    1,7,5,
    1,3,7,
    1,4,3,
    1,2,4,
    3,8,7,
    3,4,8,
    5,7,8,
    5,8,6,
    1,5,6,
    1,6,2,
    2,6,8,
    2,8,4).finished().array()-1;
}

void initSimData(sim_data& data, int num){
	cout << "initSimData\n";
	Vector3d origin = data.gridDims.cast<double>() * 0.5f;
	origin(1) *= 1.5;
	for(int i = 0; i < num; i++){
	//Vector3d offset = Vector3d::Zero();
	
		Vector3d offset = Vector3d::Random();
		offset(0) *= data.gridDims(0) * 0.25;
		offset(1) *= data.gridDims(1) * 0.25;
		offset(2) *= data.gridDims(2) * 0.25;
		data.particles.push_back(particle(origin + offset, Vector3d(0,0,0)));
		//data.particles.push_back(particle(origin + offset, Vector3d(0,1,0)+ Vector3d::Random()));
		//data.particles.push_back(particle(origin + offset, Vector3d(0,1,0) + Vector3d::Random()));
	}
	//data.particles.push_back(particle(Vector3d(-1.5, -1.5, 0), Vector3d(0.5,1,0)));
	//data.particles.push_back(particle(Vector3d(0.5, 1.5, 0), Vector3d(0.2,0.4,0)));
}

void test(){
	sim_data data(3,3,3,0.2);
	
}


int main(int argc, char *argv[])
{
  runTests();
  return 0;
  
  int gridSize = argc > 1 ? stol(argv[1]) : 3;
  double dt = argc > 2 ? stod(argv[2]) : 0.5;
  int numParts = argc > 3 ? stol(argv[3]) : 3;
  bool normalize_vel_grid = argc > 4 ? (stol(argv[4]) == 0 ? false : true) : false;
  cout << "gridSize: " << gridSize << endl;
  cout << "dt: " << dt << endl;
  cout << "numParts: " << numParts << endl;
  cout << "normalize velocity grid: " << normalize_vel_grid << endl;
  
  sim_data data(gridSize,gridSize,gridSize,0.2);
  data.normalize_velocity_grid = normalize_vel_grid;
  data.gravity = Vector3d(0,-0.4,0);
  initSimData(data, numParts);

  // Inline mesh of a cube
  double w = 1.0;
  
  MatrixXd cMin = getCubeV(w, Vector3d(0,0,0));
  MatrixXd cMax = getCubeV(w, data.gridDims.cast<double>());
  Eigen::MatrixXi f1 = getCubeF();
  Eigen::MatrixXi f2 = getCubeF() + MatrixXi::Ones(12,3)*8;
  
  Eigen::MatrixXd V = MatrixXd::Zero(16,3); //= getCubeV(w, data.gridDims.cast<double>() * 0.5);
  V << cMin, cMax;
  Eigen::MatrixXi F = MatrixXi::Zero(24,3); //= getCubeF();
  F << f1, f2;	

  // Plot the mesh
  igl::opengl::glfw::Viewer viewer;
  viewer.data().set_mesh(V, F);
  viewer.data().set_face_based(true);
  viewer.core().is_animating = true;
  
  
  //data.print_debug_once = true;
  //simulate(data, 0.017);
  viewer.callback_pre_draw = [&] (igl::opengl::glfw::Viewer&)->bool{
  		
  	// update simulation	
  	simulate(data, dt);//0.017);
    
  	// dump into matrix and set points
  	MatrixXd pts(data.particles.size(), 3);
  	for(int i = 0; i < data.particles.size(); i++){
  		pts.row(i) = data.particles[i].pos;
  	}
  	viewer.data().set_points(pts, RowVector3d(1,1,1));
  	return false;
  };
  
  viewer.callback_key_pressed = 
    [&](igl::opengl::glfw::Viewer &, unsigned int key, int mod)
  {
    double d_lambda = 0.2;
    
    Vector3d dVel = Vector3d::Zero();
    double spd = 0.5;
    if(key == 'w'){
    	dVel(1) = spd;		
    }
    else if(key == 's'){
    	dVel(1) = -spd;
    }
    else if(key == 'a'){
    	dVel(0) = -spd;
    }
    else if(key == 'd'){
    	dVel(0) = spd;
    }
    
  	data.particles[0].vel += dVel;
    switch(key)
    {
      
      case 'p':
        {
        cout << "Test key press s\n";
        data.print_debug_once = true;
        }
        break;
      case 'o':
        {
         print_data(data);
        }
        break;
      default:
        return false;
    }
    
    return true;
  };
  
  viewer.launch();
}


