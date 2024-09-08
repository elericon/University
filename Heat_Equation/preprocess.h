#ifndef PREPROCESS_H
#define PREPROCESS_H



void find_dx_dy(int nx, int ny, double &dx, double &dy){
	dx = 0.05 / nx;
	dy = 0.05 / ny;
}
void pv_imp(double *temp, int nx, int cont){
	temp[0] = cont + 1;
	temp[1] = cont + nx;
	temp[2] = cont - 1;
	temp[3] = cont - nx;
}
void mesh(int nx, int ny, double dx, double dy, double **coord, double **pv, double **bound){
	int cont = 0, cont_bound = 0;
	double alpha = 9;
	double *temp, *temp1;

	for (int j = 0; j < ny; ++j){
		for (int i = 0; i < nx; ++i){
			// find coordinates of the points
			temp = coord[cont];
			temp[0] = dx + 2 * i * dx;
			temp[1] = dy + 2 * j * dy;

			if (temp[0] <= 0.5 && temp[1] <= 0.4){
				temp[2] = 1.511e-4;
			}
			else if(temp[0] > 0.5 && temp[1] <= 0.7){
				temp[2] = 1.1364e-4;
			}
			else if(temp[0] > 0.5 && temp[1] <= 0.8){
				temp[2] = 1.2995e-4;
			}
			else if(temp[0] <= 0.5 && temp[1] <= 0.8){
				temp[2] = 6.0215e-5;
			}

			// find first neighbours
			temp1 = pv[cont];
			temp = bound[cont_bound];
			pv_imp(temp1, nx, cont);

			// no real first neighbour on the boundary
			if (i == 0){// x == 0
				temp1[2] = -1;
				temp[0] = cont;
				temp[1] = 33 * alpha / (2 * dx) ;
				if (j * dy <= 0.4){
					temp[1] /= (1500 * 750);
				}
				else{
					temp[1] /= (1900 * 810);
				}
				temp[2] = 0;//convection
				cont_bound++;
			}
			if (i == nx - 1){// x == 1.1
				temp1[0] = -1;
				temp[0] = cont;
				temp[1] = 8;
				temp[2] = 1;//pure dirichlet but in time
				cont_bound++;
			}
			if (j == 0){// y == 0
				temp1[3] = -1;
				if(i != 0 && i != nx - 1){
					temp[0] = cont;
					temp[1] = 23;
					temp[2] = 2;//pure dirichlet
					cont_bound++;
				}
			}
			if(j == ny - 1){// y == 0.8
				temp1[1] = -1;
				if(i != 0 && i != nx - 1){
					temp[0] = cont;
					temp[1] = 54.54 / (2 * dy);// dqdt * dS / dV
					if (i * dx <= 0.5){
						temp[1] /= (1900 * 810);
					}
					else{
						temp[1] /= (2500 * 930);
					}
					temp[2] = 3;//pure newman
					cont_bound++;
				}
			}
			// go to next point
			cont++;
		}
	}

	print_mat(coord,nx*ny,3,"mesh/xy.dat");
	print_mat(pv,nx*ny,4,"mesh/closest_neig.dat");
	print_mat(bound,cont_bound,3,"mesh/bound.dat");
}
void infomat_create(double **info, int nn, double **pv, double **coord, double dx, double dy){
	double *temp, *pv_i, *ci, alpha_i;
	for (int i = 0; i < nn; ++i){
		temp = info[i];
		pv_i = pv[i];
		temp[0] = 4 * dx * dy;
		alpha_i = coord[i][2];
		for (int j = 0; j < 4; ++j){
			if (pv_i[j] != -1){
				ci = coord[(int)pv_i[j]];
				if (j == 0 || j == 2){
					temp[j + 1] = 2 * dy * ((alpha_i*ci[2]) / (alpha_i + ci[2])) / dx;
				}
				else{
					temp[j + 1] = 2 * dx * ((alpha_i*ci[2]) / (alpha_i + ci[2])) / dy;
				}
			}
			else{
				temp[j + 1] = -1;
			}
		}
	}
	print_mat(info,nn,5,"mesh/info_mat.dat");
}


void preprocess(int nx, int ny){
	std::cout<<"Starting preprocess\n";
	double dx, dy;

	find_dx_dy(nx,ny,dx,dy);

	nx *= 11;
	ny *= 8;
	int nn = nx * ny;

	double **coord = mat_alloc(nn, 3);// x, y, alpha
	double **pv = mat_alloc(nn,4);// east north west south
	double **bound = mat_alloc(nn,3);//node, value, type of boundary
	
	mesh(nx,ny,dx,dy,coord,pv,bound);

	double **info = mat_alloc(nn,5), **info1;
	infomat_create(info, nn, pv, coord, dx, dy);
	std::cout<<"Preprocess completed. coord, bound, pv and info mats allocated and saved in files.\n";
}

#endif
