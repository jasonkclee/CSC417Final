#include <igl/opengl/glfw/Viewer.h>
#include <iostream>
#include "particle.h"
#include "sim_data.h"
#include "simulate.h"
#include "tests.h"
using namespace Eigen;
using namespace std;
#include <algorithm>

char* getCmdOption(char ** begin, char ** end, const std::string & option)
{
    char ** itr = std::find(begin, end, option);
    if (itr != end && ++itr != end)
    {
        return *itr;
    }
    return 0;
}

bool cmdOptionExists(char** begin, char** end, const std::string& option)
{
    return std::find(begin, end, option) != end;
}

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
	//origin(1) *= 1.5;
	for(int i = 0; i < num; i++){
	//Vector3d offset = Vector3d::Zero();
	
		Vector3d offset = Vector3d::Random();
		offset(0) *= data.gridDims(0) * 0.2;
		offset(1) *= data.gridDims(1) * 0.2;
		offset(2) *= data.gridDims(2) * 0.2;
		offset *= 0.01;
		data.particles.push_back(particle(origin + offset, Vector3d::Random()*0.001));
		//data.particles.push_back(particle(origin + offset, Vector3d(0,1,0)+ Vector3d::Random()));
		//data.particles.push_back(particle(origin + offset, Vector3d(0,1,0) + Vector3d::Random()));
	}
}




int main(int argc, char *argv[])
{
  // Variables to hold arguments
  double gravity = 0.1f;
  int gridSize = 10;
  double dt = 0.14;
  int numParts = 5;
  bool flip_method = true;
  
  if(cmdOptionExists(argv, argv+argc, "-g"))
  	gravity = stod(getCmdOption(argv, argv + argc, "-g"));
  	
  if(cmdOptionExists(argv, argv+argc, "-s"))
  	gridSize = stol(getCmdOption(argv, argv + argc, "-s"));
  
  if(cmdOptionExists(argv, argv+argc, "-dt"))
  	dt = stod(getCmdOption(argv, argv + argc, "-dt"));
  
  if(cmdOptionExists(argv, argv+argc, "-n"))
  	numParts = stol(getCmdOption(argv, argv + argc, "-n"));
  
  if(cmdOptionExists(argv, argv+argc, "-pic"))
  	flip_method = false;
  
  sim_data data(gridSize,gridSize,gridSize);
  data.gravity = Vector3d(0,-gravity,0);
  data.flip_method = flip_method;
  
  initSimData(data, numParts);

  // Cube to show boundary
  double w = 1.0;
  MatrixXd cMin = getCubeV(w, Vector3d(0,0,0));
  MatrixXd cMax = getCubeV(w, data.gridDims.cast<double>());
  Eigen::MatrixXi f1 = getCubeF();
  Eigen::MatrixXi f2 = getCubeF() + MatrixXi::Ones(12,3)*8;
  
  Eigen::MatrixXd V = MatrixXd::Zero(16,3); 
  V << cMin, cMax;
  Eigen::MatrixXi F = MatrixXi::Zero(24,3); 
  F << f1, f2;	


  igl::opengl::glfw::Viewer viewer;
  viewer.data().set_mesh(V, F);
  viewer.data().set_face_based(true);
  viewer.core().is_animating = true;
  
  
  //data.print_debug_once = true;
  viewer.callback_pre_draw = [&] (igl::opengl::glfw::Viewer&)->bool{
  		
  	// update simulation	
  	simulate(data, dt);
    
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
  	double initVel = 1.7;
    switch(key)
    {
      case 't':
      {
      	data.particles.push_back(particle(Vector3d(1.1,1.1,1.1), initVel*(Vector3d(2.0,2.0,2.0) + Vector3d::Random())));
      	break;
      }
      case 'r':
      {
        Vector3d spawn = data.gridDims.cast<double>() - Vector3d(0.1,0.1,0.1);
      	data.particles.push_back(particle(spawn, -initVel*(Vector3d(1.0,1.0,1.0) + Vector3d::Random()) ));
      	break;
      }
      case 'p':
        {
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


