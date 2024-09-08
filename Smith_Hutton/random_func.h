#ifndef RANDOM_FUNC_H
#define RANDOM_FUNC_H


double find_theta(double *temp_coord, int pv, double dx, double dy, double &S, double &v){// find velocity field at the center of the face
// return parameter to multiply the first neighbour for streamline upwind

	double prop = 1.;//parameter for the upwind part
	double x, y;

	if(pv == 0){
		x = temp_coord[0] + dx;
		y = temp_coord[1];
		S = - 2 * dy;// surface of interest
	}
	else if(pv == 1){
		x = temp_coord[0];
		y = temp_coord[1] + dy;
		S = - 2 * dx;// surface of interest
	}
	else if(pv == 2){
		x = temp_coord[0] - dx;
		y = temp_coord[1];
		S = 2 * dy;// surface of interest
	}
	else if(pv == 3){
		x = temp_coord[0];
		y = temp_coord[1] - dy;
		S = 2 * dx;// surface of interest
	}

	if(pv == 0){
		v = 2 * y * (1 - x * x);//vx
		if(v >= 0){// central node is upwind
			return 1 - prop;
		}
		else{//neighbour is upwind
			return prop;
		}
	}
	else if(pv == 1){
		v = - 2 * x * (1 - y * y);//vy
		if(v >= 0){// central node is upwind
			return 1 - prop;
		}
		else{//neighbour is upwind
			return prop;
		}
	}
	else if(pv == 2){
		v = 2 * y * (1 - x * x);//vx
		if(v >= 0){//neighbour is upwind
			return prop;
		}
		else{// central node is upwind
			return 1 - prop;
		}
	}
	else{
		v = - 2 * x * (1 - y * y);//vy
		if(v >= 0){//neighbour is upwind
			return prop;
		}
		else{// central node is upwind
			return 1 - prop;
		}
	}
}
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
