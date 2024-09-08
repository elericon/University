#ifndef PREPROCESS_H
#define PREPROCESS_H



void find_dx_dy(int nx, int ny, double &dx, double &dy){// (1/n)/2
	dx = 0.5 / nx;
	dy = 0.5 / ny;
}
void pv_imp(double *temp, int nx, int cont){
	temp[0] = cont + 1;
	temp[1] = cont + nx;
	temp[2] = cont - 1;
	temp[3] = cont - nx;
}
void mesh(int nx, int ny, double dx, double dy, double **coord, double **pv, double **bound, double RG){//GR =  gamma / rho
	int cont = 0, cont_bound = 0;
	double GR = 1 / RG;
	double alpha = 10.;
	double *temp_coord, *temp_pv, *temp_bound;

	for (int j = 0; j < ny; ++j){
		for (int i = 0; i < nx; ++i){
			// find coordinates of the points
			temp_coord = coord[cont];
			temp_coord[0] = dx + 2 * i * dx - 1;//x in [-1,1]
			temp_coord[1] = dy + 2 * j * dy;//y in [0,1]

			temp_coord[2] = GR; //"reynolds" number to check why independent of it

			// find first neighbours
			temp_pv = pv[cont];
			temp_bound = bound[cont_bound];
			pv_imp(temp_pv, nx, cont);


			// no real first neighbour on the boundary, bound[i][2] = 0,1,2,3 equates to conv, dir in time, dir, newm
			if (i == 0){// x == -1
				temp_pv[2] = -1;
				temp_bound[0] = cont;
				temp_bound[1] = 1 - tanh(alpha);
				temp_bound[2] = 2;//pure dirichlet
				temp_bound[3] = 2;// direction of the boundary of interest. 0 east, 1 north, 2 west, 3 south
				cont_bound++;
			}
			if (i == nx - 1){// x == 1 
				temp_pv[0] = -1;
				temp_bound[0] = cont;
				temp_bound[1] = 1 - tanh(alpha);
				temp_bound[2] = 2;//pure dirichlet
				temp_bound[3] = 0;// direction of the boundary of interest. 0 east, 1 north, 2 west, 3 south
				cont_bound++;
			}
			if (j == 0){// y == 0
				temp_pv[3] = -1;
				if(i <= nx / 2 && i != 0){
					temp_bound[0] = cont;
					temp_bound[1] = 1 + tanh((2 * temp_coord[0] + 1) * alpha);
					temp_bound[2] = 2;//pure dirichlet
					temp_bound[3] = 3;// direction of the boundary of interest. 0 east, 1 north, 2 west, 3 south
					cont_bound++;
				}
				else if(i != nx - 1){
					temp_bound[0] = cont;
					temp_bound[1] = 0;
					temp_bound[2] = 3;//newmann
					temp_bound[3] = 3;// direction of the boundary of interest. 0 east, 1 north, 2 west, 3 south
					cont_bound++;
				}
			}
			if (j == ny - 1){// y == 1
				temp_pv[1] = -1;
				if(i != 0 && i != nx - 1){
					temp_bound[0] = cont;
					temp_bound[1] = 1 - tanh(alpha);
					temp_bound[2] = 2;//pure dirichlet
					temp_bound[3] = 1;// direction of the boundary of interest. 0 east, 1 north, 2 west, 3 south
					cont_bound++;
				}
			}
			// go to next point
			cont++;
		}
	}

	print_mat(coord,nx*ny,3,"mesh/xy.dat");
	print_mat(pv,nx*ny,4,"mesh/closest_neig.dat");
	print_mat(bound,cont_bound,4,"mesh/bound.dat");
}
void infomat_create(double **info, int nn, double **pv, double **coord, double dx, double dy){
	double *temp, *pv_i, alpha_j, alpha_i;
	for (int i = 0; i < nn; ++i){
		temp = info[i];
		pv_i = pv[i];
		temp[0] = 4 * dx * dy;//v_i
		alpha_i = coord[i][2];
		for (int j = 0; j < 4; ++j){
			// if (pv_i[j] != -1){
			// alpha_j = coord[(int)pv_i[j]][2];
			if (j == 0 || j == 2){
				temp[j + 1] = dy * alpha_i / dx;//((alpha_i * alpha_j) / (alpha_i + alpha_j))
			}
			else{
				temp[j + 1] = dx * alpha_i / dy;//((alpha_i * alpha_j) / (alpha_i + alpha_j))
			}
			//taken away as I need the value to be stored also for the non existent boundary nodes
			// }
			// else{
			// 	temp[j + 1] = -1;
			// }
		}
		temp[5] = dx;
		temp[6] = dy;
	}
	print_mat(info,nn,7,"mesh/info_mat.dat");
}


void preprocess(int nx, int ny, double RG){
	std::cout<<"Starting preprocess\n";
	double dx, dy;

	find_dx_dy(nx,ny,dx,dy);
	nx *= 2;

	int nn = nx * ny;

	double **coord = mat_alloc(nn, 3);// x, y, alpha
	double **pv = mat_alloc(nn,4);// east north west south
	double **bound = mat_alloc(nn,4);//node, value, type of boundary
	
	mesh(nx,ny,dx,dy,coord,pv,bound,RG);

	double **info = mat_alloc(nn,7), **info1;
	infomat_create(info, nn, pv, coord, dx, dy);
	std::cout<<"Preprocess completed. coord, bound, pv and info mats allocated and saved in files.\n";
}

#endif
