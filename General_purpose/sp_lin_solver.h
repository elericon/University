#ifndef SP_LIN_SOLVER_H
#define SP_LIN_SOLVER_H


double sign(double b){
	if(b>0){
		return 1;
	}
	else if(b<0){
		return -1;
	}
	else{
		return 0;
	}
}
void triang_solver(double **M, int nr, double *b, double *x){//make sure that the matrix is upper triangolar
	double *temp;// handle for matrix
	double temp_x;//handle for x_i
	
	// Ensure M[nr-i][nr-i] is not zero before dividing
	if (M[nr-1][nr-1]!=0){
		x[nr-1] = b[nr-1]/M[nr-1][nr-1];//solve last unknown
	}

	for (int i = 2; i <= nr; ++i){
		temp = M[nr-i];
		temp_x = b[nr-i];
		for (int j = 1; j < i; ++j) {
		    temp_x -= temp[nr-j] * x[nr-j];//
		}
		// Ensure M[nr-i][nr-i] is not zero before dividing
		if (temp[nr-i]!=0){
			temp_x /= temp[nr-i]; 
		}
		else{
			std::cout << "Diagonal entry is null, problem\n";
		}
		x[nr-i] = temp_x;//called once only
	}
}

void find_sc_horiz(double *s, double *c, double *temp, int col){
	double a, b, t;
	
	a = temp[col];
	b = temp[col+1];
	
	if(abs(b)>=abs(a)){//compute the sin and cos without computational risks of overflow
		t = a/b;
		*s = sign(b)/sqrt(1+t*t);
		*c = (*s)*t;
	}
	else{
		t = b/a;
		*c = sign(a)/sqrt(1+t*t);
		*s = (*c)*t;
	}
}
// // solves sparse linear systems using the general minimization residual method
void sp_GMRES(sparse &A, double *b, double tol, long int kmax_in, double *x0, double *x, double *&rho1, int *flag, int *iter){//sparse L, sparse U,
	long int n = A.get_nrow();
	
	if (kmax_in > n){
		kmax_in = n;
		std::cout<<"Maximum number of iterations is greater than allowed by finite termination property of GMRES. Setting kmax_in equal to the size of the matrix\n";
	}

	double **v = mat_alloc(kmax_in+1,n);
	double **h = mat_alloc(kmax_in+1,kmax_in+2);//needs to be transposed at the end, use the transposed for faster column addressing
	
	double *r = new double[n];//first residual vector
	double *rho = new double[kmax_in+1];//vector of the errors
	double *temp_v = new double[n];
	double *temp_h;//handle
	double *u_j = new double[n];
	double *g = new double[n+1];
	double *sin_i = new double[n];
	double *cos_i = new double[n];

	//constants
	double beta, h_ii, nu;
	*flag = 0;//default
	
	// compute the first residual
	sparse_vm_prod(A,x0,temp_v);//temp_v=A*x
	vector_up(b,temp_v,r,n,1);//r[0] = b-temp = b - A*x0

	//compute the first norm
	beta = norm(r,n);
	if(beta + 1 == 1){
		std::cout<<"Error: beta = 0\n";
	}
	rho[0] = beta;

	// create right handside
	zero(g,n);
	g[0] = beta;

	// compute the first vector
	vector_up(r,beta,v[0],n,1);//v[0] = r/beta

	tol *= norm(b,n);

	int j;
	for (j = 0; j < kmax_in && tol < rho[j]; ++j){
		temp_h = h[j];//faster computations as one less link
		sparse_vm_prod(A,v[j],u_j);//u_j = A*vj

		for (int i = 0; i <= j; ++i){// complete the matrix row per row instead of using columns, to speed up computation
			temp_h[i] = scalar_prod(u_j,v[i],n);//h[i][j] = u_i*vj
			vector_up(v[i],temp_h[i],temp_v,n,0);//temp_v = h[i][j]*v[j]
			vector_up(u_j,temp_v,u_j,n,1);//u_i = u_i - h[i][j]*vj
		}

		h_ii = norm(u_j,n);
		temp_h[j+1] = h_ii;
	
		if(h_ii + 1 == 1){//breakdown, set flag to -1 and exit cycle
			*flag = -1;
			j++;
			std::cout<<"Lucky breakdown happened\n";
			break;
		}

		vector_up(u_j,h_ii,v[j+1],n,1);//vj+1 = u_j/h_ii
		
		for (int i = 0; i < j; ++i){
			//apply givens rotations on the trasposed h matrix
			nu = cos_i[i] * temp_h[i] + sin_i[i]* temp_h[i+1];
			temp_h[i+1] = -sin_i[i]*temp_h[i] + cos_i[i]*temp_h[i+1];
			temp_h[i] = nu;
		}

		find_sc_horiz(&sin_i[j], &cos_i[j], temp_h, j);//finds the givens factors to be used in the givens rotations

		//apply givens rotations to the transposed matrix
		temp_h[j] = cos_i[j] * temp_h[j] + sin_i[j] * temp_h[j+1];
		temp_h[j+1] = 0;

		//apply givens rotations to the RHS
		g[j+1] = -sin_i[j] * g[j];
		g[j] = cos_i[j] * g[j];

		//compute the residual
		rho[j+1] = abs1(g[j+1]); 
		if (j%100 == 0){
			std::cout<<j/kmax_in*100<<"% \n";
		}
	}
	if (j + 1 == 1){
		std::cout<<"Error: j = 0\n";
	}

	double **H = mat_alloc(j,j);
	copy_submat(H,h,0,0,j,j,1);//transpose the hessemberg matrix to solve the system
	mat_dealloc(h);//try to free the space as fast as possible before going on, as the laptop is not able to take it and kills the program

	double **V = mat_alloc(n,j);
	copy_submat(V,v,0,0,j,n,1);
	mat_dealloc(v);//try to free the space as fast as possible before going on, as the laptop is not able to take it and kills the program

	double *y = new double[j];

	triang_solver(H, j, g, y);//solve upper triangolar system and then use the found y to find the final solution

	vector_up(x0,vm_prod(V,y,n,j),x,n,0);//x = x0 + V*y
	
	rho1 = new double[j];
	for (int i = 0; i < j; ++i){
		rho1[i] = rho[i];
	}
	// convergence not reached
	if(rho[j] >= tol){
		std::cout<<"Convergence not reached: residual of "<<rho[j]<<"\n";
		*flag = 1;
	}
	else{
		std::cout<<"Convergence reached: residual of "<<rho[j]<<"\n";
	}
	
	*iter = j;


	// memory deallocation
	mat_dealloc(H);
	mat_dealloc(V);

	delete[] rho;
	delete[] r;
	delete[] temp_v;
	delete[] g;
	delete[] y;
	delete[] u_j;

	return;
}

// solves sparse linear systems using the general minimization residual method with jacobi preconditioner
void sp_GMRES(sparse &A, double *b, double tol, long int kmax_in, double *x0, double *x, double *&rho1, int *flag, int *iter, sparse M){//M preconditioner, must be jacobi
	long int n = A.get_nrow();

	if (kmax_in > n){
		kmax_in = n;
		std::cout<<"Maximum number of iterations is greater than allowed by finite termination property of GMRES. Setting kmax_in equal to the size of the matrix\n";
	}
	
	double **v = mat_alloc(kmax_in+1,n);
	double **h = mat_alloc(kmax_in+1,kmax_in+2);//needs to be transposed at the end, use the transposed for faster column addressing
	
	double *r = new double[n];//first residual vector
	double *rho = new double[kmax_in+1];//vector of the errors
	double *temp_v = new double[n];
	double *temp_h;//handle
	double *u_j = new double[n];
	double *g = new double[n+1];
	double *sin_i = new double[n];
	double *cos_i = new double[n];

	//constants
	double beta, h_ii, nu;
	*flag = 0;//default
	
	// compute the first residual
	sparse_vm_prod(A,x0,temp_v);//temp_v=A*x
	vector_up(b,temp_v,r,n,1);//r[0] = b-temp = b - A*x0

	//compute the first norm
	sparse_vm_prod(M,r,temp_v);//temp_v = M*r

	beta = norm(temp_v,n);
	if(beta + 1 == 1){
		std::cout<<"Error: beta = 0\n";
	}
	rho[0] = beta;

	// create right handside
	zero(g,n);
	g[0] = beta;

	// compute the first vector
	vector_up(temp_v,beta,v[0],n,1);//v[0] = r/beta

	double *prec_b = new double[n];
	sparse_vm_prod(M,b,prec_b);//prec_b = M*b
	tol *= norm(prec_b,n);

	int j;
	for (j = 0; j < kmax_in && tol < rho[j]; ++j){
		temp_h = h[j];//faster computations as one less link
		sparse_vm_prod(A,v[j],temp_v);//u_j = A*vj
		sparse_vm_prod(M,temp_v,u_j);//u_j = M*u_j


		for (int i = 0; i <= j; ++i){// complete the matrix row per row instead of using columns, to speed up computation
			temp_h[i] = scalar_prod(u_j,v[i],n);//h[i][j] = u_i*vj
			vector_up(v[i],temp_h[i],temp_v,n,0);//temp_v = h[i][j]*v[j]
			vector_up(u_j,temp_v,u_j,n,1);//u_i = u_i - h[i][j]*vj
		}
		
		h_ii = norm(u_j,n);
		temp_h[j+1] = h_ii;
	
		if(h_ii + 1 == 1){//breakdown, set flag to -1 and exit cycle
			*flag = -1;
			j++;
			std::cout<<"Lucky breakdown happened\n";
			break;
		}

		vector_up(u_j,h_ii,v[j+1],n,1);//vj+1 = u_j/h_ii
		
		for (int i = 0; i < j; ++i){
			//apply givens rotations on the trasposed h matrix
			nu = cos_i[i] * temp_h[i] + sin_i[i]* temp_h[i+1];
			temp_h[i+1] = -sin_i[i]*temp_h[i] + cos_i[i]*temp_h[i+1];
			temp_h[i] = nu;
		}

		find_sc_horiz(&sin_i[j], &cos_i[j], temp_h, j);//finds the givens factors to be used in the givens rotations

		//apply givens rotations to the transposed matrix
		temp_h[j] = cos_i[j] * temp_h[j] + sin_i[j] * temp_h[j+1];
		temp_h[j+1] = 0;

		//apply givens rotations to the RHS
		g[j+1] = -sin_i[j] * g[j];
		g[j] = cos_i[j] * g[j];

		//compute the residual
		rho[j+1] = abs1(g[j+1]); 
		if (j%100 == 0){
			std::cout<<(double)j/(double)kmax_in*100.<<"%, residual of "<<rho[j+1]<<" \n";
		}
	}
	if (j + 1 == 1){
		std::cout<<"Error: j = 0\n";
	}
	
	double **H = mat_alloc(j,j);
	copy_submat(H,h,0,0,j,j,1);//transpose the hessemberg matrix to solve the system
	mat_dealloc(h);//try to free the space as fast as possible before going on, as the laptop is not able to take it and kills the program

	double **V = mat_alloc(n,j);
	copy_submat(V,v,0,0,j,n,1);
	mat_dealloc(v);//try to free the space as fast as possible before going on, as the laptop is not able to take it and kills the program

	double *y = new double[j];

	triang_solver(H, j, g, y);//solve upper triangolar system and then use the found y to find the final solution

	vector_up(x0,vm_prod(V,y,n,j),x,n,0);//x = x0 + V*y
	
	rho1 = new double[j+1];
	for (int i = 0; i <= j; ++i){
		rho1[i] = rho[i];
	}

	// convergence not reached
	if(rho[j] >= tol){
		std::cout<<"Convergence not reached: residual of "<<rho[j]<<"\n";
		*flag = 1;
	}
	else{
		std::cout<<"Convergence reached in " << j+1<<" iterations. Residual of "<<rho[j]<<"\n";
	}

	*iter = j + 1;

	// memory deallocation
	mat_dealloc(H);
	mat_dealloc(V);
	

	delete[] rho;
	delete[] r;
	delete[] temp_v;
	delete[] g;
	delete[] y;
	delete[] u_j;

	return;
}
void lower_triang_solver(double **M, int nr, double *b, double *x){//make sure that the matrix is lower triangolar
	double *temp;// handle for matrix
	double temp_x;//handle for x_i
	
	// Ensure M[0][0] is not zero before dividing
	if (M[0][0]!=0){
		x[0] = b[0]/M[0][0];//solve last unknown
	}

	for (int i = 1; i < nr; ++i){
		temp = M[i];
		temp_x = b[i];
		for (int j = 0; j < i; ++j) {
		    temp_x -= temp[j] * x[j];//
		}
		// Ensure M[nr-i][nr-i] is not zero before dividing
		if (temp[i]!=0){
			temp_x /= temp[i]; 
		}
		else{
			std::cout << "Diagonal entry is null, problem\n";
		}
		x[i] = temp_x;//called once only
	}
}
void LU_decomposition(double **M, double **L, double **U, int n) {
    //use handles for speed
    double *temp_U, *temp_L, *temp_M;

    // Initialize L to identity and suppose already allocated to zero
    #pragma omp parallel for num_threads(n_threads)
    for (int i = 0; i < n; i++) {
            L[i][i] = 1.0;
    }

    #pragma omp parallel for num_threads(n_threads)
    for (int k = 0; k < n; k++) {
    	temp_U = U[k];
    	temp_M = M[k];
    	temp_L = L[k];

        // Upper triangular matrix U
        for (int j = k; j < n; j++) {
            temp_U[j] = temp_M[j];

            for (int p = 0; p < k; p++) {
                temp_U[j] -= temp_L[p] * U[p][j];
            }
        }
        // Lower triangular matrix L
        for (int i = k + 1; i < n; i++) {
        	temp_L = L[i];

            temp_L[k] = M[i][k];
            for (int p = 0; p < k; p++) {
                temp_L[k] -= temp_L[p] * U[p][k];
            }
            temp_L[k] /= temp_U[k];
        }
    }
}
void solve_LU(double **L, double **U, double *rhs, long int nn_p, double *res){
	double *partial_res = new double[nn_p]();

	lower_triang_solver(L, nn_p, rhs, partial_res);
	triang_solver(U, nn_p, partial_res, res);

	delete[] partial_res;
}

void keep_last(double **V, int nr, int nc){
	double *V_last = new double[nc]();

	copy_vector(V_last,V[nr - 1],nc);
	zero(V,nr,nc);
	copy_vector(V[0],V_last,nc);

	delete[] V_last;
}

// solves sparse linear systems using the general minimization residual method with jacobi preconditioner
void restarted_GMRES(sparse &A, double *b, double tol, long int kmax_in, double *x0, double *x, double *&rho1, int &flag, int &iter, sparse M){//M preconditioner, must be jacobi
	long int n = A.get_nrow();

	int rest_in = 50;

	if (kmax_in > n){
		kmax_in = n;
		std::cout<<"Maximum number of iterations is greater than allowed by finite termination property of GMRES. Setting kmax_in equal to the size of the matrix\n";
	}
	
	double **v = mat_alloc(rest_in+1,n);
	double **h = mat_alloc(rest_in+1,rest_in+2);//needs to be transposed at the end, use the transposed for faster column addressing
	
	double *r = new double[n];//first residual vector
	double *rho = new double[3 * kmax_in / rest_in + 1];//vector of the errors
	double *temp_v = new double[n];
	double *temp_h;//handle
	double *u_j = new double[n];
	double *g = new double[rest_in+1];
	double *sin_i = new double[n];
	double *cos_i = new double[n];

	//constants
	double beta, h_ii, nu;
	flag = 0;//default
	iter = 0;

	// compute the first residual
	sparse_vm_prod(A,x0,temp_v);//temp_v=A*x
	vector_up(b,temp_v,r,n,1);//r[0] = b-temp = b - A*x0

	//compute the first norm
	sparse_vm_prod(M,r,temp_v);//temp_v = M*r

	beta = norm(temp_v,n);
	if(beta + 1 == 1){
		std::cout<<"Error: beta = 0\n";
	}
	rho[0] = beta;
	// create right handside
	zero(g,n);
	g[0] = beta;

	// compute the first vector
	vector_up(temp_v,beta,v[0],n,1);//v[0] = r/beta

	double *prec_b = new double[n];
	sparse_vm_prod(M,b,prec_b);//prec_b = M*b
	tol *= norm(prec_b,n);

	int j_in;
	for (int j_out = 0; j_out < 3 * kmax_in / rest_in && tol < rho[iter]; ++j_out){//outer loop

		if (j_out != 0){//keep last and zero the rest
			keep_last(v,rest_in+1,n);
			zero(h,rest_in+1,rest_in+2);
			zero(g,n+1);
			g[0] = rho[iter];
		}

		for (j_in = 0; j_in < rest_in && tol < rho[j_in + iter]; ++j_in){//inner loop
			temp_h = h[j_in];//faster computations as one less link
			sparse_vm_prod(A,v[j_in],temp_v);//u_j = A*vj
			sparse_vm_prod(M,temp_v,u_j);//u_j = M*u_j
	
	
			for (int i = 0; i <= j_in; ++i){// complete the matrix row per row instead of using columns, to speed up computation
				temp_h[i] = scalar_prod(u_j,v[i],n);//h[i][j] = u_i*vj
				vector_up(v[i],temp_h[i],temp_v,n,0);//temp_v = h[i][j]*v[j]
				vector_up(u_j,temp_v,u_j,n,1);//u_i = u_i - h[i][j]*vj
			}
			

			h_ii = norm(u_j,n);
			temp_h[j_in+1] = h_ii;
		
			if(h_ii + 1 == 1){//breakdown, set flag to -1 and exit cycle
				flag = -1;
				j_in++;
				std::cout<<"Lucky breakdown happened\n";
				break;
			}
	
			vector_up(u_j,h_ii,v[j_in+1],n,1);//vj+1 = u_j/h_ii
			
			for (int i = 0; i < j_in; ++i){
				//apply givens rotations on the trasposed h matrix
				nu = cos_i[i] * temp_h[i] + sin_i[i]* temp_h[i+1];
				temp_h[i+1] = -sin_i[i]*temp_h[i] + cos_i[i]*temp_h[i+1];
				temp_h[i] = nu;
			}
	
			find_sc_horiz(&sin_i[j_in], &cos_i[j_in], temp_h, j_in);//finds the givens factors to be used in the givens rotations
	
			//apply givens rotations to the transposed matrix
			temp_h[j_in] = cos_i[j_in] * temp_h[j_in] + sin_i[j_in] * temp_h[j_in+1];
			temp_h[j_in+1] = 0;
	
			//apply givens rotations to the RHS
			g[j_in+1] = -sin_i[j_in] * g[j_in];
			g[j_in] = cos_i[j_in] * g[j_in];
	
			//compute the residual
			rho[j_in + iter + 1] = abs1(g[j_in+1]); 
		}

		iter += j_in;
	}
	if (iter + 1 == 1){
		std::cout<<"Error: j = 0\n";
	}
	double **H = mat_alloc(j_in,j_in);
	copy_submat(H,h,0,0,j_in,j_in,1);//transpose the hessemberg matrix to solve the system
	double **V = mat_alloc(n,j_in);
	copy_submat(V,v,0,0,j_in,n,1);
	
	double *y = new double[j_in];
	triang_solver(H, j_in, g, y);//solve upper triangolar system and then use the found y to find the final solution

	vector_up(x0,vm_prod(V,y,n,j_in),x,n,0);//x = x0 + V*y
	
	rho1 = new double[iter+1];
	for (int i = 0; i <= iter; ++i){
		rho1[i] = rho[i];
	}

	// convergence not reached
	if(rho[iter] >= tol){
		std::cout<<"Convergence not reached: residual of "<<rho[iter]<<"\n";
		flag = 1;
	}
	else{
		std::cout<<"Convergence reached in " << iter+1<<" iterations. Residual of "<<rho[iter]<<"\n";
	}
	iter++;

	// memory deallocation
	std::cout<<"test\n";
	mat_dealloc(H);
	std::cout<<"test\n";
	mat_dealloc(V);
	std::cout<<"test\n";
	mat_dealloc(v);
	std::cout<<"test\n";
	mat_dealloc(h);
	std::cout<<"test\n";

	delete[] rho;
	delete[] r;
	delete[] temp_v;
	delete[] g;
	delete[] y;
	delete[] u_j;

	return;
}
// solves sparse linear systems using the general minimization residual method with jacobi preconditioner
void pcg(sparse &A, double *b, double tol, double *x0, double *x, double *&rho1, int &iter){
	long int n = A.get_nrow();

	copy_vector(x,x0,n);
	
	double *r = new double[n]();
	double *rho = new double[n+1]();//vector of the errors
	double *res = new double[n]();
	double *temp_v = new double[n]();
	double *z = new double[n]();
	double *p = new double[n]();
	double alpha, beta;
	
	// compute the first residual
	sparse_vm_prod(A,x0,temp_v);//temp_v=A*x
	vector_up(b,temp_v,r,n,1);//r = b-temp = b - A*x0
	copy_vector(p,r,n);//p = r

	rho[0] = scalar_prod(r,p,n);
	res[0] = norm(r,n);

	tol *= norm(b,n);//preconditioned norm

	int j;
	for (j = 0; j < n && res[j] >= tol; ++j){

		//compute the first norm
		sparse_vm_prod(A,p,z);//z = A*p
		// std::cout<<"Z\n";
		// print_vett(z,n);

		alpha = rho[j] / scalar_prod(z,p,n);//norm_r / ((A*p)^T*p)
		// std::cout<<"Alpha = "<<alpha<<"\n";

		// std::cout<<"X\n";
		lin_comb(x,1.,p,alpha,x,n);//x = x + alpha * p
		// print_vett(x,n);

		// std::cout<<"R\n";
		lin_comb(r,1.,z,-alpha,r,n);//r = r - alpha * A * p
		// print_vett(r,n);

		rho[j + 1] = scalar_prod(r,r,n);
		beta = rho[j + 1] / rho[j];
		// std::cout<<"Beta = "<<beta<<"\n";

		// std::cout<<"P\n";
		lin_comb(r,1.,p,beta,p, n);//p = r + beta * p
		// print_vett(p,n);

		res[j + 1] = sqrt(rho[j + 1]);
	}

	iter = j + 1;
	if (rho[j] < tol){
		// std::cout << "Convergence reached in " << iter <<" iterations. Residual of " << res[j] << "\n";
	}
	else{
		std::cout << "Convergence not reached: residual of " << res[j] << "\n";
	}

	rho1 = new double[iter];
	for (int i = 0; i < iter; ++i){
		rho1[i] = res[i];
	}
	

	return;
}

// solves sparse linear systems using the general minimization residual method with jacobi preconditioner
void pcg(sparse &A, double *b, double tol, double *x0, double *x, double *&rho1, int &iter, sparse M){//M preconditioner, must be jacobi
	long int n = A.get_nrow();

	copy_vector(x,x0,n);
	
	double *r = new double[n]();
	double *rho = new double[n+1]();//vector of the errors
	double *res = new double[n]();
	double *temp_v = new double[n]();
	double *z = new double[n]();
	double *p = new double[n]();
	double alpha, beta;
	
	// compute the first residual
	sparse_vm_prod(A,x0,temp_v);//temp_v=A*x
	vector_up(b,temp_v,r,n,1);//r = b-temp = b - A*x0
	sparse_vm_prod(M,r,p);//p = M * r

	rho[0] = scalar_prod(r,p,n);
	res[0] = norm(r,n);

	tol *= norm(b,n);//preconditioned norm

	int j;
	for (j = 0; j < n && res[j] >= tol; ++j){

		//compute the first norm
		sparse_vm_prod(A,p,z);//z = A*p

		alpha = rho[j] / scalar_prod(z,p,n);//norm_r / ((A*p)^T*p)
		lin_comb(x,1.,p,alpha,x,n);//x = x + alpha * p
		lin_comb(r,1.,z,-alpha,r,n);//r = r - alpha * A * p

		sparse_vm_prod(M,r,temp_v);//temp_v = M * r
		
		rho[j + 1] = scalar_prod(r,temp_v,n);
		beta = rho[j + 1] / rho[j];
		
		lin_comb(temp_v,1.,p,beta,p, n);//p = M * r + beta * p
		res[j + 1] = norm(r,n);
	}

	iter = j + 1;
	if (rho[j] < tol){
		// std::cout << "Convergence reached in " << iter <<" iterations. Residual of " << res[j] << "\n";
	}
	else{
		std::cout << "Convergence not reached: residual of " << res[j] << "\n";
	}

	rho1 = new double[iter];
	for (int i = 0; i < iter; ++i){
		rho1[i] = res[i];
	}

	delete[] r;
	delete[] rho;
	delete[] res;
	delete[] temp_v;
	delete[] z;
	delete[] p;

	return;
}

// solves sparse linear systems using the general minimization residual method with jacobi preconditioner
void pcg(sparse &A, double *b, double tol, double *x0, double *x, int &iter, sparse M){//M preconditioner, must be jacobi
	long int n = A.get_nrow();

	copy_vector(x,x0,n);
	
	double *r = new double[n]();
	double *rho = new double[n+1]();//vector of the errors
	double *res = new double[n]();
	double *temp_v = new double[n]();
	double *z = new double[n]();
	double *p = new double[n]();
	double alpha, beta;
	
	// compute the first residual
	sparse_vm_prod(A,x0,temp_v);//temp_v=A*x
	vector_up(b,temp_v,r,n,1);//r = b-temp = b - A*x0
	sparse_vm_prod(M,r,p);//p = M * r

	rho[0] = scalar_prod(r,p,n);
	res[0] = norm(r,n);

	tol *= norm(b,n);//preconditioned norm

	int j;
	for (j = 0; j < n && res[j] >= tol; ++j){

		//compute the first norm
		sparse_vm_prod(A,p,z);//z = A*p

		alpha = rho[j] / scalar_prod(z,p,n);//norm_r / ((A*p)^T*p)
		lin_comb(x,1.,p,alpha,x,n);//x = x + alpha * p
		lin_comb(r,1.,z,-alpha,r,n);//r = r - alpha * A * p

		sparse_vm_prod(M,r,temp_v);//temp_v = M * r
		
		rho[j + 1] = scalar_prod(r,temp_v,n);
		beta = rho[j + 1] / rho[j];
		
		lin_comb(temp_v,1.,p,beta,p, n);//p = M * r + beta * p
		res[j + 1] = norm(r,n);
	}

	iter = j + 1;
	if (rho[j] < tol){
		// std::cout << "Convergence reached in " << iter <<" iterations. Residual of " << res[j] << "\n";
	}
	else{
		std::cout << "Convergence not reached: residual of " << res[j] << "\n";
	}

	// rho1 = new double[iter];
	// for (int i = 0; i < iter; ++i){
	// 	rho1[i] = res[i];
	// }

	delete[] r;
	delete[] rho;
	delete[] res;
	delete[] temp_v;
	delete[] z;
	delete[] p;

	return;
}



#endif
