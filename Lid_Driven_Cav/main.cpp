#include <iostream>
#include <fstream>
#include <cmath>
#include <random>

int numcall = 0;
const int n_threads = 2;
double U_ref = 1;
double L = 1;
int unif = 1;//0 for uniform, 1 for non uniform
int randomm = 0;//0 quasi zero initial condition, 1 random initial condition

#include "omp.h"
#include "mat_vec_funct.h"
#include "mat_allocation.h"
#include "sparse_mat.h"
#include "poiss_div.h"
#include "sp_lin_solver.h"
#include "non_unif_preprocess.h"
#include "conv_diff.h"

using namespace std;

void initial_cond(double *u, double *v, long int nn_u, long int nn_v);
void print_for_postprocess(long int nn_u, long int nn_v, long int nn_p, int k, int save_cont, double **u_save, double **v_save, double **p_save, double *u, double *v, double *p);


int main(){
	double Re = 10000;
	int nx = 100, ny = 100;//nx, ny number of points in a line of 1 m
	double rho = 1;
	double mu = rho*L*U_ref/Re;
	int time = 1;// one to save time dependent solution

	//construct mesh and info matrices, as well as boundary conditions and first neighbours
	double deltat = preprocess(nx,ny,mu,rho);
	// cout<<"For next tries use deltat = "<<deltat<<" or less\n";
	// deltat /= 2;

	long int nn_u, nb_u, nn_v, nb_v, nn_p, nc;
	double final_time = 1000;
	long int M = final_time / deltat;
	int save_cont = 0, cont_iter;

	//import matrices from preprocessing
	double **bound_u = import_mat("mesh/bound_u.dat",nb_u, nc);
	double **pv_u = import_mat("mesh/pv_u.dat", nn_u, nc);
	double **info_u = import_mat("mesh/info_mat_u.dat", nn_u, nc);

	double **bound_v = import_mat("mesh/bound_v.dat",nb_v, nc);
	double **pv_v = import_mat("mesh/pv_v.dat", nn_v, nc);
	double **info_v = import_mat("mesh/info_mat_v.dat", nn_v, nc);

	double **pv_P = import_mat("mesh/pv_p.dat", nn_p, nc);
	double **info_P = import_mat("mesh/info_mat_p.dat", nn_p, nc);

	//allocate space for the solution, one for the current time step and one for the previous
	double **u = mat_alloc(2,nn_u), *temp_u;
	double **v = mat_alloc(2,nn_v), *temp_v;
	double **p = mat_alloc(2,nn_p), *temp_p;

	//save only 500 time steps out of all the result. too much space otherwise
	int n_save = 500;
	double **u_save = mat_alloc(n_save,nn_u);
	double **v_save = mat_alloc(n_save,nn_v);
	double **p_save = mat_alloc(n_save,nn_p);

	//as in article to check oscillations for Re = 10000
	long int monitor_node_u = (long int)((nx + 1) * ny * 13. / 16. + (nx + 1) * 14. / 16.);
	long int monitor_node_v = (long int)(nx * ny * 13. / 16. + nx * 14. / 16.);
	double **monitoring_point = mat_alloc(M, 2);

	cout<<"matrices allocated\n";

	//allocate space for R_u, R_v for 2 iterations
	double **R_u = mat_alloc(2, nn_u);
	double **R_v = mat_alloc(2, nn_v);

	//allocate space for the trial velocities and the poisson rhs
	double *u_trial = new double[nn_u](), *R_uEff = new double[nn_u]();
	double *v_trial = new double[nn_v](), *R_vEff = new double[nn_v]();
	double *poiss_rhs = new double[nn_p]();
	
	cout<<"variables allocated\n";

	
	// //impose initial condition, done as suppose all zero initial condition
	initial_cond(u[0],v[0],nn_u, nn_v);
	
	//build the poisson sparse matrix one as it is always the same
	sparse Poiss;
	poisson_mat(Poiss,nn_p,pv_P);
	cout<<"Computed Poisson matrix\n";
	

	//prepare for pcg
	double *residual, *u_res = new double[nn_u]();
	int iter, k = M - 1; 

	//time solution with expl euler
	for (int i = 0; i < M - 1; ++i){//M - 1
		if(norm(u_res,nn_u)>2e6){
			cout<<"diverging\n";
			k = i + 1;
			break;
		}

		if (i%100 == 0){
			cout<<"Starting time iteration number = " << i+1<<"/"<<M<<"\n";
		}

		//compute R_u, R_v using previous time step result
		R_uv(rho,nx,ny,u[0],R_u[1],pv_u,info_u,bound_u,nn_u,nb_u,
					   v[0],R_v[1],pv_v,info_v,bound_v,nn_v,nb_v);
		
		//compute the adams bashford term for R
		if(i > 0){
			compute_R_uv_eff(R_u,R_uEff,R_v,R_vEff,nn_u,nn_v);
		}
		else{
			copy_vector(R_uEff,R_u[1],nn_u);
			copy_vector(R_vEff,R_v[1],nn_v);
		}
		
		//compute trial velocities
		compute_trial_vel(deltat, rho, nx, info_u, nn_u, u_trial, u[0], R_uEff,
                                           info_v, nn_v, v_trial, v[0], R_vEff);
		
		//solve poisson equation
		poiss_RHS(rho,deltat,nn_p,info_P,pv_P,u_trial,v_trial,poiss_rhs);
				
		// cout << "Starting Poisson solution of step " << i+1 << "\n";
		pcg(Poiss, poiss_rhs, 1e-9, p[0], p[1], iter, Poiss.diag_jac());

		//compute the new velocities
		new_vel(deltat/rho, p[1], nx, ny, info_u, u_trial, u[1], 
										  info_v, v_trial, v[1]);
		
		if (i%100 == 0){
			cout<<"u_res norm = "<<norm(u_res,nn_u)<<", mean iteration last 100 = "<<cont_iter/100<<"\n";
			cont_iter = 0;
		}
		else{
			cont_iter += iter;
		}
		// check_div_free(pv_P,u[1],v[1],info_P,nn_p,i);
		monitoring_point[i][0] = u[1][monitor_node_u];
		monitoring_point[i][1] = v[1][monitor_node_v];

		//check convergence
		vector_up(u[1],u[0],u_res,nn_u,1);
		if ((norm(u_res,nn_u)<5e-8|| iter<5) && i > 20){
			k = i + 1;
			cout<<"u_res norm = "<<norm(u_res,nn_u)<<", i = "<<i<<", iter = "<<iter<<"\n";
			break;
		}
		
		if (i % ((M/10)/n_save + 1) == 0 && save_cont < n_save){
			copy_vector(u_save[save_cont],u[1],nn_u);
			copy_vector(v_save[save_cont],v[1],nn_v);
			copy_vector(p_save[save_cont],p[1],nn_p);
			save_cont++;
		}
		else if(save_cont > n_save){
			cout<<"saved too much, not until the end.\n";
		}

		//update R for the AB iteration
		copy_vector(R_u[0], R_u[1], nn_u);
		copy_vector(R_v[0], R_v[1], nn_v);
		copy_vector(u[0],u[1],nn_u);
		copy_vector(v[0],v[1],nn_v);
		copy_vector(p[0],p[1],nn_p);

	}
	
	//free space
	mat_dealloc(pv_u);
	mat_dealloc(pv_v);
	mat_dealloc(pv_P);

	// //free space
	mat_dealloc(bound_u);
	mat_dealloc(bound_v);
	
	//free space
	mat_dealloc(info_u);
	mat_dealloc(info_v);
	mat_dealloc(info_P);

	print_mat(monitoring_point,k,2,"monitor.dat");
	mat_dealloc(monitoring_point);

	print_for_postprocess(nn_u,nn_v,nn_p,k,save_cont,u_save,v_save,p_save,u[1],v[1],p[1]);

	mat_dealloc(u_save);
	mat_dealloc(v_save);
	mat_dealloc(p_save);

	mat_dealloc(u);
	mat_dealloc(v);
	mat_dealloc(p);

	mat_dealloc(R_u);
	mat_dealloc(R_v);

	delete[] u_trial;
	delete[] v_trial;
	delete[] poiss_rhs;
	delete[] R_vEff;
	delete[] R_uEff;
	delete[] u_res;


	return 0;
}



void initial_cond(double *u, double *v, long int nn_u, long int nn_v){
	if (randomm == 1){
		for (int i = 0; i < nn_u; ++i){
			u[i] = random_11();
		}
		for (int i = 0; i < nn_v; ++i){
			v[i] = random_11();
		}
	}
	else{//not exactly zero otherwise problems
		for (int i = 0; i < nn_v; ++i){
			v[i] = random_11() / 1e17;//
		}
		// for (int i = 0; i < nn_u; ++i){
		// 	u[i] = 1;//random_11() / 1e17
		// }
	}
}

void print_for_postprocess(long int nn_u, long int nn_v, long int nn_p, int k, int save_cont, double **u_save, double **v_save, double **p_save, double *u, double *v, double *p){
	long int nn = max(nn_u, nn_v), nc;
	ofstream result1, result2, result3, result4;
	int nt = min(k, save_cont);
	double *temp_u, *temp_v, *temp_p;

	
	//print results on a file to be animated on matlab
	result1.open("Results/u_in_time.dat");
	result2.open("Results/v_in_time.dat");
	result3.open("Results/p_in_time.dat");
	
	for (int i = 0; i < save_cont; i++){
		temp_u = u_save[i];
		temp_v = v_save[i];
		temp_p = p_save[i];
		
		for (int j = 0; j < nn; ++j){
			if (j < nn_u){
				result1 << temp_u[j] << "\t";
			}
			if (j < nn_v){
				result2 << temp_v[j] << "\t";
			}
			if (j < nn_p){
				result3 << temp_p[j] << "\t";
				}
		}

		result1<<"\n";
		result2<<"\n";
		result3<<"\n";	
	}

	result1.close();
	result2.close();
	result3.close();
	
	//print results on a file to be animated on matlab
	result1.open("Results/u_final.dat");
	result2.open("Results/v_final.dat");
	result3.open("Results/p_final.dat");

	for (int j = 0; j < nn; ++j){
		if (j < nn_u){
			result1 << u[j] << "\t";
		}
		if (j < nn_v){
			result2 << v[j] << "\t";
		}
		if (j < nn_p){
			result3 << p[j] << "\t";
		}
	}
	result1<<"\n";
	result2<<"\n";
	result3<<"\n";	

	result1.close();
	result2.close();
	result3.close();

	double **coord_u = import_mat("mesh/xy_u.dat",nn_u, nc);
	double **coord_v = import_mat("mesh/xy_v.dat",nn_v, nc);
	double **coord_P = import_mat("mesh/xy_p.dat",nn_p, nc);
	
	//get results on the mid line
	result1.open("Results/u_mid.dat");
	result2.open("Results/v_mid.dat");

	for (int i = 0; i < nn; ++i){
		if (i < nn_u){
			// if(coord_u[i][0] == 0.5){
			if(coord_u[i][0] <= 0.5 + 0.005 && coord_u[i][0] >= 0.5 - 0.005){ 
				result1 << coord_u[i][1] << "\t" << u[i] << "\n";
			}
		}
		if (i < nn_v){
			// if(coord_v[i][1] == 0.5){
			if(coord_v[i][1] <= 0.5 + 0.005 && coord_v[i][1] >= 0.5 - 0.005){
				result2 << coord_v[i][0] << "\t" << v[i] << "\n";
			}
		}
	}
	
	result1.close();
	result2.close();

	// double *v_mod = vel_mod(u, v, nn_p);
	// print_vett(v_mod,nn_p,"Results/v_mod.dat");


	//free space
	mat_dealloc(coord_u);
	mat_dealloc(coord_v);
	mat_dealloc(coord_P);
}
