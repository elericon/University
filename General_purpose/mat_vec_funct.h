#ifndef MAT_VEC_FUNCT_H
#define MAT_VEC_FUNCT_H

// #include "omp.h"
extern int numcall;

double norm(double *v, int n){//computes the modulus of a given vector
	double m = 0., t;

	#pragma omp parallel for reduction(+:m) num_threads(n_threads)
    for (int i = 0; i < n; ++i){
        t = v[i];
		m += t*t;
    }
 
	return sqrt(m);
}
double max(double *D, int n){
	double max = D[0];
    for (int i = 0; i < n; ++i){
        if(D[i] > max){
			max = D[i];
		}
    }
	return max;
}
double max(double x, double y){
	if(x >= y){
		return x;
	}
	else{
		return y;
	}
}
double min(double x, double y){
	if(x <= y){
		return x;
	}
	else{
		return y;
	}
}
double min(double *D, int n){
	double min = D[0];
    for (int i = 0; i < n; ++i){
        if(D[i] < min){
			min = D[i];
		}
    }
	return min;
}
double sum(double *v, int n){
	double sum = 0.;

	for (int i = 0; i < n; ++i){
		sum += v[i];
	}
	return sum;
}
double abs1(double n){
	if (n>=0){
		return n;
	}
	else{
		return (-n);
	}
}
void mean_zero(double *v, int nn){
	double mean = 0;

	#pragma omp parallel for reduction(+:mean) num_threads(n_threads)
	for (int i = 0; i < nn; ++i){
		mean += v[i];
	}
	mean /= nn;

	#pragma omp parallel for num_threads(n_threads)
	for (int i = 0; i < nn; ++i){
		v[i] -= mean;
	}
}
double mean(double *v, int nn){
	double mean = 0;
	#pragma omp parallel for reduction(+:mean) num_threads(n_threads)
	for (int i = 0; i < nn; ++i){
		mean += v[i];
	}
	mean /= nn;

	return mean;
}
double abs_mean(double *v, int nn){
	double mean = 0;
	for (int i = 0; i < nn; ++i){
		mean += abs1(v[i]);
	}
	mean /= nn;

	return mean;
}
void print_mat(double **M, int n){//prints a full matrix
    for (int i = 0; i < n; ++i){
        for (int j = 0; j < n; ++j){
			std::cout<<M[i][j]<< "\t";
		}
		std::cout<<std::endl;
    }
	std::cout<<std::endl;
}
void print_mat(double **M, int n, std::string filename){//prints matrix on file
	std::ofstream print;
	print.open(filename);
	for (int i = 0; i < n; ++i){
        for (int j = 0; j < n; ++j){
			print<<M[i][j]<< "\t";
		}
		print<<std::endl;
    }
	print<<std::endl;
	print.close();
}

void print_mat(double **M, int nr, int nc){//prints a full matrix
    for (int i = 0; i < nr; ++i){
        for (int j = 0; j < nc; ++j){
			std::cout<<M[i][j]<< "\t";
		}
		std::cout<<std::endl;
    }
	std::cout<<std::endl;
}
void print_mat(double **M, int nr, int nc, std::string filename){//prints matrix on file
	std::ofstream print;
	print.open(filename);
	
	print << nr << "\t" << nc;
	for (int i = 0; i < nc - 2; ++i){
		print << "\t" << 0;
	}
	print << "\n";
	for (int i = 0; i < nr; ++i){
        for (int j = 0; j < nc; ++j){
			print<<M[i][j]<< "\t";
		}
		print<<std::endl;
    }
	print<<std::endl;
	print.close();
}
void print_vett(double *v, int n, std::string filename){//prints vector on file
	std::ofstream print;
	print.open(filename);
	
	for (int i = 0; i < n; ++i){
        print<<v[i]<< "\t";
    }
	print<<std::endl;
	print.close();
}
void print_vett(double *M, int n){//prints a full vector
    for (int i = 0; i < n; ++i){
        std::cout<<M[i]<< "\t";
    }
	std::cout<<"\n\n";
}
void print_mat(int **M, int nr, int nc){//prints a full matrix
    for (int i = 0; i < nr; ++i){
        for (int j = 0; j < nc; ++j){
			std::cout<<M[i][j]<< "\t";
		}
		std::cout<<std::endl;
    }
	std::cout<<std::endl;
}
void print_vett(int *M, int n){//prints a full vector
    for (int i = 0; i < n; ++i){
        std::cout<<M[i]<< "\t";
    }
	std::cout<<"\n\n";
}

void copy_vector(int *res, int *v, int n){
	#pragma omp parallel for num_threads(n_threads)
    for (int i = 0; i < n; ++i){
      	res[i] = v[i];
    }
}
void copy_vector(double *res, double *v, int n){
	#pragma omp parallel for num_threads(n_threads)
    for (int i = 0; i < n; ++i){
      	res[i] = v[i];
    }
}
void zero(double *v, int n){
	#pragma omp parallel for num_threads(n_threads)
    for (int i = 0; i < n; ++i){
       	v[i] = 0;
    }
}
void zero(int *v, int n){
	#pragma omp parallel for num_threads(n_threads)
    for (int i = 0; i < n; ++i){
        v[i] = 0;
    }
}
void zero(double **M, int nr, int nc) {
	double *temp;
	#pragma omp parallel for private(temp) num_threads(n_threads)
    for (int i = 0; i < nr; ++i){   
    	temp = M[i];  
        for (int j = 0; j < nc; ++j) {   
         	temp[j] = 0;
        }
    }
}
void transpose(double **A, double **At, int nr, int nc){
	double *temp;
	#pragma omp parallel for private(temp) num_threads(n_threads)
	for (int i = 0; i < nr; ++i){
		temp = A[i];
		for (int j = 0; j < nc; ++j){
			At[j][i] = temp[j];
		}
	}
}
void copy_submat(double **Result,double **A,int nr_r_start, int nc_r_start,int nr_r_end, int nc_r_end, int opt){//opt = 0 copy the submat, opt = 1 copy the transposed submat
	double *temp;
	if (opt == 0){

		#pragma omp parallel for private(temp) num_threads(n_threads)
		for (int i = nr_r_start; i < nr_r_end; ++i){
			temp = A[i];
			for (int j = nc_r_start; j < nc_r_end; ++j){
				Result[i-nr_r_start][j-nc_r_start] = temp[j];
			}
		}
	}
	else{

		#pragma omp parallel for private(temp) num_threads(n_threads)
		for (int i = nr_r_start; i < nr_r_end; ++i){
			temp = A[i];
			for (int j = nc_r_start; j < nc_r_end; ++j){
				Result[j-nc_r_start][i-nr_r_start] = temp[j];
			}
		}
	}
}
void vector_up(double *x, double *y, double *res, int n, bool type){// updates the vector in setting res as the result of the operation
	if(type == 0){
		#pragma omp parallel for num_threads(n_threads)
		for (int i = 0; i < n; ++i){
			res[i] = x[i] + y[i];
		}
	}
	else{
		#pragma omp parallel for num_threads(n_threads)
		for (int i = 0; i < n; ++i){
			res[i] = x[i] - y[i];
		}
	}
	return;
}
void lin_comb(double *x, double alpha, double *y, double beta, double *res, int n){
	#pragma omp parallel for num_threads(n_threads)
	for (int i = 0; i < n; ++i){
		res[i] = alpha * x[i] + beta * y[i];
	}
}
void vector_up(double *x, double y, double *res, int n, bool type){// updates the vector in setting res as the result of the operation
	if(type == 0){

		#pragma omp parallel for num_threads(n_threads)
		for (int i = 0; i < n; ++i){
			res[i] = x[i]*y;
		}
	}
	else{

		#pragma omp parallel for num_threads(n_threads)
		for (int i = 0; i < n; ++i){
			res[i] = x[i]/y;
		}
	}
	return;
}
void vector_up(int *x, int *y, int *res, int n, bool type){// updates the vector in setting res as the result of the operation
	if(type == 0){

		#pragma omp parallel for num_threads(n_threads)
		for (int i = 0; i < n; ++i){
			res[i] = x[i] + y[i];
		}
	}
	else{

		#pragma omp parallel for num_threads(n_threads)
		for (int i = 0; i < n; ++i){
			res[i] = x[i] - y[i];
		}
	}
	return;
}

int nnz(double **M, int row){
	// find the number of nonzeros or else said n_term
	int cont = 0;
	double *temp;//handle
	// #pragmma omp parallel for private(temp)
	for (int i = 0; i < row; ++i){
		temp = M[i];
		for (int j = 0; j < row; ++j){
			if(temp[j] + 1 != 1){//exclude terms under the machine epsilon
		
				// #pragmma omp atomic
				cont++;
			}
		}
	}
	return cont;
}
int nnz(double **M, int row, int col){
	// find the number of nonzeros or else said n_term
	int cont = 0;
	double *temp;//handle
	// #pragmma omp parallel for private(temp)
	for (int i = 0; i < row; ++i){
		temp = M[i];
		for (int j = 0; j < col; ++j){
			if(temp[j] + 1 != 1){//exclude terms under the machine epsilon
		
				// #pragmma omp atomic
				cont++;
			}
		}
	}
	return cont;
}
int nnz_vec(int *M, int n){
	// find the number of nonzeros or else said n_term
	int cont = 0;
	// #pragmma omp parallel for
	for (int i = 0; i < n; ++i){
		if(M[i] + 1 != 1){//exclude terms under the machine epsilon
	
			// #pragmma omp atomic
			cont++;
		}
	}
	return cont;
}

double *vm_prod(double **A, double *x, int nr, int nc){
	double *b = new double[nr];
	double *temp;

	#pragma omp parallel for private(temp) num_threads(n_threads)
	for (int i = 0; i < nr; ++i){
		b[i] = 0;
		temp = A[i];
		for (int j = 0; j < nc; ++j){
			b[i] += temp[j]*x[j];
		}
	}
	return b;
}
double scalar_prod(double *v1, double *v2, int n){
	double s = 0.;
	
	#pragma omp parallel for reduction(+:s) num_threads(n_threads)
	for (int i = 0; i < n; ++i){
		s += v1[i]*v2[i];
	}
	return s;
}
void create_ej(int n, int j, double *w){
	zero(w,n);
	w[j] = 1;
}

void ones(double *v, int n){
	#pragma omp parallel for num_threads(n_threads)
	for (int i = 0; i < n; ++i){
		v[i] = 1;
	}
}

void vec_assemble(double *res, double *x, double *y, int nx, int ny){
	for (int i = 0; i < nx; ++i){
		res[i] = x[i];
	}
	for (int i = 0; i < ny; ++i){
		res[i + nx] = y[i];
	}
}

#endif
