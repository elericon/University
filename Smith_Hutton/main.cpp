#include <iostream>
#include <fstream>
#include <cmath>

int numcall = 0;
const int n_threads = 1;

#include "mat_vec_funct.h"
#include "mat_allocation.h"
#include "sparse_mat.h"
#include "random_func.h"
#include "sp_lin_solver.h"
#include "preprocess.h"
#include "conv_diff.h"

using namespace std;

int main(){
	int nx = 50, ny = 50;//nx, ny number of points in a line of 1 m
	double RG = 1e6;
	cout<<"RG = "<<RG<<"\n";
	int time = 1;// zero for time independent and 1 for time dependent

	//construct mesh and info matrices, as well as boundary conditions and first neighbours
	preprocess(nx,ny,RG);
	nx *= 2;

	int nn, nc, nb;
	int final_time = 3;
	double deltat = 0.0001;
	long int M = final_time / deltat;

	//import matrices from preprocessing
	double **coord = import_mat("mesh/xy.dat",nn, nc);
	double **bound = import_mat("mesh/bound.dat",nb, nc);
	double **pv = import_mat("mesh/closest_neig.dat", nn, nc);
	double **info = import_mat("mesh/info_mat.dat", nn, nc);

	cout<<"matrices allocated\n";


	//allocate space for the solution
	double **phi = mat_alloc(M,nn), *temp = new double[nn]();

	cout<<"phi allocated\n";

	
	//impose initial condition
	double *phi0 = phi[0];
	for (int i = 0; i < nn; ++i){
		phi0[i] = 0;
	}

	double *resv = new double[nn];
	int k;

	if(time == 0){
		sparse H, C;
		double *q = new double[nn](), *rho = new double[nn]();
		int flag = 0, iter;

		create_CD_mat(H, C, coord, pv, info, nn);

		H.sumv(C.get_coef());

		bound_cond(bound, nb, H, q);
		cout<<"starting GMRES \n";
		sp_GMRES(H, q, 1e-7, 20*sqrt(8*nn), phi[0], phi[M - 1], rho, &flag, &iter, H.diag_jac());
	}
	else if(time == 1){
		//time solution
		for (int i = 0; i < M - 1; ++i){
			create_CD_RHS(phi[i], temp, coord, pv, info, bound, nn, nb);//rhs(temp) = (-conv+diff)/(rho*dV)
			vector_up(temp, deltat, temp, nn, 0);//temp = rhs * dt
			vector_up(temp, phi[i], phi[i + 1], nn, 0);// update phi[i+1]
		}
	}
	

	if( time == 1){//norm(resv, nn) < 1000 &&

		//print results on a file to be animated on matlab
		ofstream result1;
		result1.open("result_in_time.dat");
		for (int i = 0; i < M; ++i){
			for (int j = 0; j < nn; ++j){
				result1<<phi[i][j]<<"\t";
			}
			result1<<"\n";	
		}
		result1.close();
	}


	//print results on a file to be animated on matlab
	ofstream result;
	result.open("result_outlet.dat");
	phi0 = phi[M - 1];

	for (int i = 0; i < nx; ++i){
		result << coord[i][0] << "\t" << coord[i][1] << "\t" << phi0[i] << "\n";
	}
	result.close();
	

	//print results on a file to be animated on matlab
	ofstream result2;
	result2.open("result_final.dat");
	phi0 = phi[M - 1];

	for (int i = 0; i < nn; ++i){
		result2 << phi0[i] << "\t";
	}
	result2 << "\n";
	result2.close();
	

	

	return 0;
}
