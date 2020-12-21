#include "sim_data.h"
void build_grid(std::vector<Eigen::MatrixXd>& grid, int w, int h, int l){
	for(int z = 0; z < l; z++){
		grid.push_back(Eigen::MatrixXd(w,h));
	}
}
