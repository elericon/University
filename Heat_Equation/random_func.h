#ifndef RANDOM_FUNC_H
#define RANDOM_FUNC_H



// void dir_bound(double **coord, int nn){
// 	double tol = 1e-8;
// 	double *temp, *temp1, x, y;
// 	int nb;

// 	std::ofstream dir;
// 	dir.open("mesh/bound_dir005v1.dat");

// 	dir >> nb;
// 	double **bound = mat_alloc(nb,2);

// 	for (int j = 0; j < nb; ++j){
// 		temp = bound[j];
// 		bnd = (int)temp[0]-1;//boundary node
		
// 		temp1 = coord[bnd];
// 		x = temp[0];
// 		y = temp[1];

// 		if (abs1(x - 0.)<tol){
// 			count++;
// 			dir<<i+1<<"\n";
// 		}
// 		else if(abs1(x - 1.1)<tol){
// 			count++;
// 			dir<<i+1<<"\n";
// 		}
// 		else if(abs1(y - 0.)<tol){
// 			count++;
// 			dir<<i+1<<"\n";
// 		}
// 		else if(abs1(y - 0.8)<tol){
// 			count++;
// 			dir<<i+1<<"\n";
// 		}
// 		else{}

// 	}
// 	std::cout<<"Number of boundary nodes is Nb = "<<count<<"\n";
// 	dir.close();
// }


// boundary condition imposition with penalty method
void bound_cond(double **bound, long int nb, sparse &H, double *q){
	double R = 1e15;
	int *iat = H.get_iat();
	int *Ja = H.get_Ja();
	double *v_coef = H.get_coef();
	double *temp;
	int bnd, k, iat0, iat1; 

	#pragma omp parallel for private(temp, iat0, iat1, bnd, k) num_threads(n_threads)
	for (int j = 0; j < nb; ++j){

		temp = bound[j];
		bnd = (int)temp[0];//boundary node

		if (temp[2] != 3 && temp[2] != 0){
			iat0 = iat[bnd];
			iat1 = iat[bnd+1];
			
			for (k = iat0; k < iat1; ++k){
				if (Ja[k] == bnd){
					break;
				}
			} 
			v_coef[k] = R;
			q[bnd] = R * (temp[1]);//R*boundary value (dirichlet and dirichlet in time)
		}
		else if(temp[2] == 0){
			iat0 = iat[bnd];
			iat1 = iat[bnd+1];
			
			for (k = iat0; k < iat1; ++k){
				if (Ja[k] == bnd && temp[2] == 0){//impose convective effects on the matrix
					v_coef[k] += temp[1] / 33;
					break;
				}
			} 
			q[bnd] += temp[1];//impose convective rhs
		}
		else{
			q[bnd] = temp[1];//impose newmann rhs
		}
	}
}
// boundary condition imposition with penalty method
// as H1, H2 have the same distribution I just need to modify the value they are equal to
void bound_cond(double **bound, long int nb, sparse &H1, sparse &H2, double *q){
	double R = 1e15;
	int *iat = H1.get_iat();
	int *Ja = H1.get_Ja();
	int nn = H1.get_nrow();
	double *v_coef1 = H1.get_coef(), *v_coef2 = H2.get_coef();
	double *temp;
	int bnd, k, iat0, iat1; 

	#pragma omp parallel for private(temp, iat0, iat1, bnd, k) num_threads(n_threads)
	for (int j = 0; j < nb; ++j){

		temp = bound[j];
		bnd = (int)temp[0]-1;//boundary node
		
		iat0 = iat[bnd];
		iat1 = iat[bnd+1];
		
		for (k = iat0; k < iat1; ++k){
			if (Ja[k] == bnd){
				break;
			}
		} 
		v_coef1[k] = R;
		v_coef2[k] = R;
		q[bnd] = R*(temp[1]);//R*boundary value
		q[bnd + nn] = R*(temp[2]);//R*boundary value
	}
}

#endif
