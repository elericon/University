#ifndef CONV_DIFF_H
#define CONV_DIFF_H


void create_CD_mat(sparse &H_sparse, sparse &C_sparse, double **coord, double **pv, double **info, long int nn){// implicit
    double **H = mat_alloc(nn,nn), *tempH, *tempC, *tempI, *pvi, temp, Vi, S, theta, v;
    double **C = mat_alloc(nn,nn), *temp_coord, vx, vy, dx = info[0][5], dy = info[0][6];
    int pv_int;

    //full matrix assembly
    for (int i = 0; i < nn; ++i){
    	tempH = H[i];
        tempC = C[i];
    	tempI = info[i];
    	pvi = pv[i];
        temp_coord = coord[i];

    	for (int j = 0; j < 4; ++j){
    		Vi = tempI[0];// finite volume of element i
            pv_int = (int)pvi[j];

    		if (pv_int != -1){
                // diffusion part
    			temp = tempI[j + 1] / Vi;
    			tempH[pv_int] = temp;
    			tempH[i] -= temp;

                // convection part implicit
                // compute the velocity field at the given finite volume interface and surface of interest
                theta = find_theta(temp_coord, j, dx, dy, S, v);// how much upwind to keep
                tempC[pv_int] = theta * v * S / Vi;
                tempC[i] -= (1 - theta) * v * S / Vi;
    		}
    	}
    }

	//render the matrix sparse for both H and C
	H_sparse.alloc(nn,nnz(H,nn));
	render_sparse(H, &H_sparse);
    mat_dealloc(H);

    C_sparse.alloc(nn,nnz(C,nn));
    render_sparse(C, &C_sparse);
    mat_dealloc(C);
}
void create_CD_RHS(double *res_n, double *rhs, double **coord, double **pv, double **info, double **bound, long int nn, int nb){
//res_n is the result at iteration n of the algorithm, rhs needs to be zeroed
    double *tempI, *pvi, *temp_coord, temp, Vi, theta, S, v, res_ni, *bnd;
    double vx, vy, dx = info[0][5], dy = info[0][6], phi_n, D;
    double *diff = new double[nn](), *conv = new double[nn]();
    int pv_int, bnd_node;

    for (int i = 0; i < nn; ++i){
        tempI = info[i];// information matrix of point i
        pvi = pv[i];// first neighbours of point i
        temp_coord = coord[i];// coord of point i
        res_ni = res_n[i];// result of time n at point i

        Vi = tempI[0];// finite volume of element i

        for (int j = 0; j < 4; ++j){
            pv_int = (int)pvi[j];

            if (pv_int != -1){
                // diffusion part
                diff[i] += tempI[j + 1] * (res_n[pv_int] - res_ni) / Vi;

                // convection part
                // compute the velocity field at the given finite volume interface and surface of interest
                theta = find_theta(temp_coord, j, dx, dy, S, v);// how much upwind to keep. at the moment 1/2 for the upwind and 1/2 for the downwind
                conv[i] += ((1 - theta) * res_ni + theta * res_n[pv_int]) * v * S / Vi;   

            }
        }
    }
    //boundary conditions
    for (int k = 0; k < nb; ++k){
        bnd = bound[k];
        bnd_node = (int)bnd[0];
        temp_coord = coord[bnd_node];// coord of point bnd_node
        tempI = info[bnd_node];// information matrix of point bnd_node
        res_ni = res_n[bnd_node];// result of time n at point bnd_node
        Vi = tempI[0];// finite volume of element bnd_node
        int j = (int)bnd[3];

        if (bnd[2] == 2){//pure dirichlet non time dependent
            diff[bnd_node] += 2 * tempI[j + 1] * (bnd[1] - res_ni) / Vi;

            // compute the velocity field at the given finite volume interface and surface of interest
            theta = find_theta(temp_coord, j, dx, dy, S, v);// how much upwind to keep. at the moment 1/2 for the upwind and 1/2 for the downwind
            conv[bnd_node] += bnd[1] * v * S / Vi;   
        }
        else if(bnd[2] == 3){// pure neumann
            if (bnd[3] == 0 || bnd[3] == 2){
                phi_n = bnd[1] * dy + res_ni; 
                D = dy;
            }
            else{
                phi_n = bnd[1] * dx + res_ni; 
                D = dx;
            }

            diff[bnd_node] += 2 * D * tempI[j + 1] * bnd[1] / Vi;

            // compute the velocity field at the given finite volume interface and surface of interest
            theta = find_theta(temp_coord,j,dx,dy,S,v);// how much upwind to keep. at the moment 2/3 for the upwind and 1/3 for the downwind
            conv[bnd_node] += res_ni * v * S / Vi;  
        }
    }
    vector_up(diff, conv, rhs, nn, 0); //conv+diff

    // Clean up dynamically allocated memory
    delete[] diff;
    delete[] conv;
}
#endif
