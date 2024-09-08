#ifndef POISS_DIV_H
#define POISS_DIV_H


double div_pointP(double *pvi, double *u, double *v, double Sx, double Sy){//returns the divergence computed at point P 
    int pv_int;
    double div = 0.;

    //cycle over first neighbours
    for (int j = 0; j < 4; ++j){
        pv_int = (int)pvi[j];

        if (j == 0){//east
            div += u[pv_int] * Sx;
        }
        else if (j == 1){//north
            div += v[pv_int] * Sy;
        }
        else if (j == 2){//west
            div -= u[pv_int] * Sx;
        }
        else{//south
            div -= v[pv_int] * Sy;
        }
    }
    return div;
}
void poisson_mat(sparse &Poiss_sparse, long int nn_p, double **pv_p){
    
    double *pvi, *po_i;
    double **Poiss = mat_alloc(nn_p,nn_p);
    int pv_int;

    #pragma omp parallel for num_threads(n_threads)
    for (int i = 0; i < nn_p; ++i){
        pvi = pv_p[i];
        po_i = Poiss[i];

        for (int j = 0; j < 4; ++j){
            pv_int = (int)pvi[j + 4];

            if (pv_int != -1){
                po_i[pv_int] = -1;
                po_i[i]++; 
            }
        }
    }
    
    Poiss[nn_p - (int)sqrt(nn_p) - 1][nn_p - (int)sqrt(nn_p) - 1] = 1e8;//set the pressure in point nn_p-1 to zero

    //render the matrix sparse
    Poiss_sparse.alloc(nn_p,nnz(Poiss,nn_p));
    render_sparse(Poiss, &Poiss_sparse);
    mat_dealloc(Poiss);
}

void poiss_RHS(double rho, double dt, long int nn_p, double **info_p, double **pv_P, double *u_trial, double *v_trial, double *rhs){//rhs already allocated
    // zero(rhs, nn_p);//to be sure for no errors

    double *pvi, rhot = - rho / dt, *tempIP;
    // std::cout<<"Kx = "<<rhot*Sx<<", Ky = "<<rhot*Sy<<"\n";

    //cycle over all the nodes
    #pragma omp parallel for num_threads(n_threads)
    for (int i = 0; i < nn_p; ++i){
        pvi = pv_P[i];
        tempIP = info_p[i];
        rhs[i] = rhot * div_pointP(pvi,u_trial,v_trial,tempIP[6],tempIP[5]);
    }
}

void check_div_free(double **pv_p, double *u_n1, double *v_n1, double **info_p, long int nn_p, int iter){
    double *pvi, Sx, Sy, *res = new double[nn_p](), *tempIP;
    
    //cycle over all the nodes
    #pragma omp parallel for num_threads(n_threads)
    for (int i = 0; i < nn_p; ++i){
        tempIP = info_p[i];
        pvi = pv_p[i];
        res[i] = div_pointP(pvi,u_n1,v_n1,tempIP[6],tempIP[5]);
    }
    if(iter%100 == 0){
        // print_vett(res,nn_p,"Results/div.dat");
        std::cout<<"divergence norm = "<<norm(res,nn_p)<<"\n";
        std::cout<<"divergence mean = "<<mean(res,nn_p)<<"\n";
        std::cout<<"divergence abs_mean = "<<abs_mean(res,nn_p)<<"\n";
        std::cout<<"\n";
    }

    delete[] res;
}

extern double U_ref;

// Function to generate a random number between -1 and 1
double random_11() {
    static std::random_device rd; // Seed the random number generator
    static std::mt19937 gen(rd()); // Mersenne Twister engine
    static std::uniform_real_distribution<> dis(-1.0, 1.0); // Define the range

    return dis(gen); // Generate and return a random number
}



#endif
