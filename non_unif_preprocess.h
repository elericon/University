#ifndef NON_UNIF_PREPROCESS_H
#define NON_UNIF_PREPROCESS_H

extern double U_ref;
extern double L;
extern int unif;

void find_dx_dy(int &nx, int &ny, double *dx, double *dy) {
    if(nx%2 == 1){//so I have a node on the half line to draw the velocity easily
        nx++;
        std::cout << "Modified nx to have a node in the center.\n";
    }
    if(ny%2 == 1){//so I have a node on the half line to draw the velocity easily
        ny++;
        std::cout << "Modified ny to have a node in the center.\n";
    }
    // Define total length in x and y directions
    double Lx = L;
    double Ly = L;

    // Define a stretching factor
    double stretch_factor;
    if (unif == 0){
        stretch_factor = 1.;
    }
    else{
        stretch_factor = 2.;
    }
    
    // Calculate total sum for normalization
    double sum_x = 0.0, sum_y = 0.0;
    for (int i = 0; i < nx; ++i) {
        double xi = (double)i / (nx - 1); // Normalized position
        dx[i] = std::pow(stretch_factor, 1.0 - abs1(1.0 - 2.0 * xi));
        sum_x += dx[i];
    }
    for (int j = 0; j < ny; ++j) {
        double yj = (double)j / (ny - 1); // Normalized position
        dy[j] = std::pow(stretch_factor, 1.0 - abs1(1.0 - 2.0 * yj));
        sum_y += dy[j];
    }
    //uso formula prof, poi traslo di meta distanza x0 e x0 diventa il primo dx, poi iterativamente trovo i successivi. idem in y

    // Normalize to ensure sum(dx) = Lx and sum(dy) = Ly
    for (int i = 0; i < nx; ++i) {
        dx[i] *= Lx / sum_x;
    }
    for (int j = 0; j < ny; ++j) {
        dy[j] *= Ly / sum_y;
    }
}



void pv_imp(double *temp, int i, int j, int nx, int ny, int uvp){
    //impose those from the different meshes
    if (uvp == 0){//u: top right then anticlock, all from v
        temp[0] = i + (nx - 1) * (j + 1);
        temp[1] = i - 1 + (nx - 1) * (j + 1);
        temp[2] = i - 1 + (nx - 1) * j;
        temp[3] = i + (nx - 1) * j;

        if (i == 0){//left side
            temp[1] = -1;
            temp[2] = -1;
        }
        if (i == (nx - 1)){//right side
            temp[0] = -1;
            temp[3] = -1;
        }
    }
    else if (uvp == 1){//v: top right then anticlock, all from u
        temp[0] = i + 1 + (nx + 1) * j;
        temp[1] = i + (nx + 1) * j;
        temp[2] = i + (nx + 1) * (j - 1);
        temp[3] = i + 1 + (nx + 1) * (j - 1);

        if (j == ny){//top side
            temp[0] = -1;
            temp[1] = -1;
        }
        if (j == 0){//bottom side
            temp[2] = -1;
            temp[3] = -1;
        }
    }
    else{//p: E,N,W,S, 0 and 2 from u, 1 and 3 from v
        temp[0] = i + 1 + (nx + 1) * j;//from u
        temp[1] = i + nx * (j + 1);//from v
        temp[2] = i + (nx + 1) * j;//from u
        temp[3] = i + nx * j;//from v
    }

    //impose those from the same mesh
    double cont = i + nx * j;
    temp[4] = cont + 1;
    temp[5] = cont + nx;
    temp[6] = cont - 1;
    temp[7] = cont - nx;
}
void imp_BC(double *temp_bound, int cont, double val, int type, int dir, int &cont_bound){
    temp_bound[0] = cont;
    temp_bound[1] = val;//value of the constrained point
    temp_bound[2] = type;//pure dirichlet
    temp_bound[3] = dir;// direction of the boundary of interest. 0 east, 1 north, 2 west, 3 south
    cont_bound++;
}
void tri_mesh(int nx, int ny, double *dx, double *dy, double mu) {
    long int nn_p = nx * ny, nn_u = (nx + 1) * ny, nn_v = nx * (ny + 1);
    double **pv_u = mat_alloc(nn_u, 8);// the first 4 from the other meshes, last 4 from the same mesh
    double **coord_u = mat_alloc(nn_u, 3);
    double **bound_u = mat_alloc(nn_u, 4);//node, value, type of boundary, direction of boundary

    double **pv_v = mat_alloc(nn_v, 8);// the first 4 from the other meshes, last 4 from the same mesh
    double **coord_v = mat_alloc(nn_v, 3);
    double **bound_v = mat_alloc(nn_v, 4);//node, value, type of boundary, direction of boundary

    double **pv_p = mat_alloc(nn_p, 8);// the first 4 from the other meshes, last 4 from the same mesh
    double **coord_p = mat_alloc(nn_p, 3);

    int cont_u = 0, cont_v = 0, cont_P = 0;
    int cont_bound_u = 0, cont_bound_v = 0; 

    double *temp_coord_u, *temp_pv_u, *temp_bound_u;
    double *temp_coord_v, *temp_pv_v, *temp_bound_v;
    double *temp_coord_p, *temp_pv_p;

    // Initialize cumulative coordinates
    double x_cumulative = 0.0;
    double y_cumulative = 0.0;

    // nx++;
    // ny++;
    
    //create mesh and boundary conditions for all the variables
    for (int j = 0; j < ny + 1; ++j){
        x_cumulative = 0.0; // Reset x cumulative for each row
        for (int i = 0; i < nx + 1; ++i){
            // find coordinates of the points staggered mesh for u,v,P
            if (j < ny){
                temp_coord_u = coord_u[cont_u];
                temp_coord_u[0] = x_cumulative;//x
                temp_coord_u[1] = y_cumulative + dy[j] / 2.0;//y
                temp_coord_u[2] = mu; //"reynolds" number 
                
                // find first neighbours for u
                temp_pv_u = pv_u[cont_u];
                temp_bound_u = bound_u[cont_bound_u];
                pv_imp(temp_pv_u,i,j,nx+1,ny,0);//the +1 is to correct the number of points in x
            }
            if (i < nx){
                temp_coord_v = coord_v[cont_v];
                // std::cout<<"dx_v = "<<dx[i]<<", dy_v = "<<dy[j]<<"\n";
                temp_coord_v[0] = x_cumulative + dx[i] / 2.0;//x
                temp_coord_v[1] = y_cumulative;//y
                temp_coord_v[2] = mu; //"reynolds" number 

                // find first neighbours for v
                temp_pv_v = pv_v[cont_v];
                temp_bound_v = bound_v[cont_bound_v];
                pv_imp(temp_pv_v,i,j,nx,ny,1);//the +1 is to correct the number of points in y
            }
            if (i < nx && j < ny){
                temp_coord_p = coord_p[cont_P];
                // std::cout<<"dx_p = "<<dx[i]<<", dy_p = "<<dy[j]<<"\n";
                temp_coord_p[0] = x_cumulative + dx[i] / 2.0;//x
                temp_coord_p[1] = y_cumulative + dy[j] / 2.0;//y
                temp_coord_p[2] = mu; //"reynolds" number 

                // find first neighbours for P
                temp_pv_p = pv_p[cont_P];
                pv_imp(temp_pv_p,i,j,nx,ny,2);
            }

            // no real first neighbour on the boundary, bound[i][2] = 0,1,2,3 equates to conv, dir in time, dir, newm
            if (j < ny){//u
                if (i == 0){// left side
                    temp_pv_u[6] = -1;//from same mesh
                    if(j != ny - 1 && j != 0){//otherwise conflicting/repeated boundary conditions
                        imp_BC(temp_bound_u, cont_u, 0, 2, 2, cont_bound_u);
                    }
                }
                if (i == nx){// right side, boundary cond only for u
                    temp_pv_u[4] = -1;//from same mesh
                    if(j != ny - 1 && j != 0){//otherwise conflicting/repeated boundary conditions
                        imp_BC(temp_bound_u, cont_u, 0, 2, 0, cont_bound_u);
                    }
                }
                if (j == 0){// bottom side
                    temp_pv_u[7] = -1;//from same mesh
                    imp_BC(temp_bound_u, cont_u, 0, 2, 3, cont_bound_u);
                }
                if (j == ny - 1){// top side
                    temp_pv_u[5] = -1;//from same mesh
                    imp_BC(temp_bound_u, cont_u, U_ref, 2, 1, cont_bound_u);
                }

                // go to next point
                cont_u++;
            }
            if (i < nx){//v
                if (i == 0){// left side
                    temp_pv_v[6] = -1;//from same mesh
                    if(j != ny && j != 0){//otherwise conflicting/repeated boundary conditions
                        imp_BC(temp_bound_v, cont_v, 0, 3, 2, cont_bound_v);
                    }
                }
                if (i == nx - 1){// right side, boundary cond only for v
                    temp_pv_v[4] = -1;//from same mesh
                    if(j != ny && j != 0){//otherwise conflicting/repeated boundary conditions
                        imp_BC(temp_bound_v, cont_v, 0, 3, 0, cont_bound_v);
                    }
                }
                if (j == 0){// bottom side
                    temp_pv_v[7] = -1;//from same mesh
                    imp_BC(temp_bound_v, cont_v, 0, 3, 3, cont_bound_v);
                }
                if (j == ny){// top side
                    temp_pv_v[5] = -1;//from same mesh
                    imp_BC(temp_bound_v, cont_v, 0, 3, 1, cont_bound_v);
                }

                // go to next point
                cont_v++;
            }
            if (i < nx && j < ny){//P
                // no real first neighbour on the boundary, bound[i][2] = 0,1,2,3 equates to conv, dir in time, dir, newm
                if (i == 0){// left side
                    temp_pv_p[6] = -1;//from same mesh
                }
                if (i == nx - 1){// right side
                    temp_pv_p[4] = -1;//from same mesh
                }
                if (j == 0){// bottom side
                    temp_pv_p[7] = -1;//from same mesh
                }
                if (j == ny - 1){// top side
                    temp_pv_p[5] = -1;//from same mesh
                }
                // go to next point
                cont_P++;
            }
            // Update cumulative coordinates
            if (i < nx+1) x_cumulative += dx[i];

        }
        if (j < ny+1) y_cumulative += dy[j];
    }


    // print matrices on external files to split in preprocess and process
    print_mat(coord_u,nn_u,3,"mesh/xy_u.dat");
    print_mat(pv_u,nn_u,8,"mesh/pv_u.dat");
    print_mat(bound_u,cont_bound_u,4,"mesh/bound_u.dat");

    print_mat(coord_v,nn_v,3,"mesh/xy_v.dat");
    print_mat(pv_v,nn_v,8,"mesh/pv_v.dat");
    print_mat(bound_v,cont_bound_v,4,"mesh/bound_v.dat");

    print_mat(coord_p,nn_p,3,"mesh/xy_p.dat");
    print_mat(pv_p,nn_p,8,"mesh/pv_p.dat");


    //free matrices
    mat_dealloc(coord_u);
    mat_dealloc(pv_u);
    mat_dealloc(bound_u);

    mat_dealloc(coord_v);
    mat_dealloc(pv_v);
    mat_dealloc(bound_v);

    mat_dealloc(coord_p);
    mat_dealloc(pv_p);
}

void diff_factor(double *pvi, double *tempI, double lambda, double *dx, double *dy, int colx, int rowy, int nx, int ny){
    double dist;
    int pv_int;

    for (int j = 0; j < 4; ++j) {
        pv_int = (int)pvi[j + 4];
        if (j == 0 || j == 2) {
            if(pv_int != -1){
                dist = (dx[colx % nx] + dx[pv_int % nx]) / 2;
            }
            else{
                dist = dx[colx % nx] / 2;
            }
            tempI[j + 1] = dy[rowy % ny] * lambda / dist;
        } 
        else {
            if(pv_int != -1){
                dist = (dy[rowy % ny] + dy[pv_int]) / 2;
            }
            else{
                dist = dy[rowy % ny] / 2;
            }
            tempI[j + 1] = dx[colx % nx] * lambda / dist;
        }
    }
}
void infomat_create(int nx, int ny, double *dx, double *dy) {
    long int nn_P, nn_u, nn_v, nc;

    double **pv_u = import_mat("mesh/pv_u.dat", nn_u, nc);
    double **pv_v = import_mat("mesh/pv_v.dat", nn_v, nc);
    double **pv_p = import_mat("mesh/pv_p.dat", nn_P, nc);
    
    double **coord_u = import_mat("mesh/xy_u.dat", nn_u, nc);
    double **coord_v = import_mat("mesh/xy_v.dat", nn_v, nc);
    double **coord_P = import_mat("mesh/xy_p.dat", nn_P, nc);

    double **info_P = mat_alloc(nn_P, 7), **info_u = mat_alloc(nn_u, 7), **info_v = mat_alloc(nn_v, 7);

    double *tempIU, alphaU_i, *pvi_u;
    double *tempIV, alphaV_i, *pvi_v;
    double *tempIP, alphaP_i, *pvi_p;

    for (int j = 0; j < ny + 1; ++j) {
        for (int i = 0; i < nx + 1; ++i) {
            //p
            if (i < nx && j < ny) { //check not out of bounds
                pvi_p = pv_p[i + nx * j];
                tempIP = info_P[i + nx * j];
                tempIP[0] = dx[i] * dy[j]; //V_i
                alphaP_i = coord_P[i + nx * j][2];
                tempIP[5] = dx[i];
                tempIP[6] = dy[j];
                diff_factor(pvi_p, tempIP, alphaP_i, dx, dy, i, j, nx, ny);
            }
            //u
            if (j < ny) { //check not out of bounds
                pvi_u = pv_u[i + (nx + 1) * j];
                tempIU = info_u[i + (nx + 1) * j];
                alphaU_i = coord_u[i + (nx + 1) * j][2];
                tempIU[6] = dy[j];
                if (i != nx){
                    tempIU[0] = dx[i] * dy[j]; //V_i
                    tempIU[5] = dx[i];
                    diff_factor(pvi_u, tempIU, alphaU_i, dx, dy, i, j, nx+1, ny);
                }
                else{
                    tempIU[0] = dx[i-1] * dy[j]; //V_i
                    tempIU[5] = dx[i-1];
                    diff_factor(pvi_u, tempIU, alphaU_i, dx, dy, i-1, j, nx+1, ny);
                }
                
            }
            //v
            if (i < nx) {
                pvi_v = pv_v[i + nx * j];
                tempIV = info_v[i + nx * j];
                alphaV_i = coord_v[i + nx * j][2];
                tempIV[5] = dx[i];
                if (j != ny){
                    tempIV[0] = dx[i] * dy[j]; //V_i
                    tempIV[6] = dy[j];
                    diff_factor(pvi_v, tempIV, alphaV_i, dx, dy, i, j, nx, ny+1);
                }
                else{
                    tempIV[0] = dx[i] * dy[j-1]; //V_i
                    tempIV[6] = dy[j-1];
                    diff_factor(pvi_v, tempIV, alphaV_i, dx, dy, i, j-1, nx, ny+1);
                }
            }
        }
    }

    print_mat(info_u, nn_u, 7, "mesh/info_mat_u.dat");
    print_mat(info_v, nn_v, 7, "mesh/info_mat_v.dat");
    print_mat(info_P, nn_P, 7, "mesh/info_mat_p.dat");

    mat_dealloc(coord_u);
    mat_dealloc(coord_v);
    mat_dealloc(coord_P);

    mat_dealloc(pv_u);
    mat_dealloc(pv_v);
    mat_dealloc(pv_p);

    mat_dealloc(info_u);
    mat_dealloc(info_v);
    mat_dealloc(info_P);
}



double preprocess(int nx, int ny, double mu, double rho) {
    std::cout << "Starting preprocess\n";
    double *dx = new double[nx](), *dy = new double[ny]();
    double dt_conv, dt_diff, delta;

    find_dx_dy(nx, ny, dx, dy);
    tri_mesh(nx, ny, dx, dy, mu);
    infomat_create(nx, ny, dx, dy);

    // Compute deltat to satisfy the Courant conditions
    delta = min(min(dx,nx), min(dy,ny));
    dt_conv = 0.3 * delta / U_ref;
    dt_diff = 0.15 * rho * delta * delta / mu;

    std::cout << "Preprocess completed. coord, bound, pv and info mats allocated and saved in files for u,v and p.\n";
    return min(dt_conv, dt_diff);
}



#endif