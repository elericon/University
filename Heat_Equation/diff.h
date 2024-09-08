#ifndef DIFF_H
#define DIFF_H


void create_D_mat(sparse &H_sparse, double **coord, double **pv, double **info, long int nn){
    double **H = mat_alloc(nn,nn), *tempH, *tempI, *pvi, temp, Vi;

    //full matrix assembly
    for (int i = 0; i < nn; ++i){
    	tempH = H[i];
    	tempI = info[i];
    	pvi = pv[i];

    	for (int j = 0; j < 4; ++j){
    		Vi = tempI[0];// finite volume of element i

    		if (pvi[j] != -1){// convection part
    			temp = tempI[j + 1] / Vi;
    			tempH[(int)pvi[j]] = temp;
    			tempH[i] -= temp;
    		}
    	}
    }

	//render the matrix sparse
	H_sparse.alloc(nn,nnz(H,nn));
	render_sparse(H, &H_sparse);
    mat_dealloc(H);
}
#endif
