#ifndef CONV_DIFF_H
#define CONV_DIFF_H

extern double U_ref;

double diff_uv(double *uv_n, double uv_P, double *tempI, double *pvi){
    double diff = 0;
    int pv_int;

    for (int j = 4; j < 8; ++j){//take the first neighbours from u or v, depending on which it is passed
        pv_int = (int)pvi[j];

        if (pv_int != -1){
            // diffusion part
            diff += tempI[j - 3] * (uv_n[pv_int] - uv_P);
        }
    }

    return diff;
}
double conv_uu(double *u_n, double *v_n, double *pv_u, double rho, double Sx, double Sy, double u_P){//use CDS convective scheme
    double *u = new double[6]();// u_e, u_n, u_w, u_s, v_n, v_s
    double *u_prov = new double[8]();
    double check = 1e15;
    int pv_int;

    for (int i = 0; i < 4; ++i){//find velocities from v
        pv_int = (int)pv_u[i];
        if (pv_int != -1){
            u_prov[i] = v_n[pv_int];
        }
        else{
            u_prov[i] = check;//big enough to see big errors if not handled correctly
        }
    }
    for (int i = 4; i < 8; ++i){//find velocities from u
        pv_int = (int)pv_u[i];
        if (pv_int != -1){
            u_prov[i] = u_n[pv_int];
        }
        else{
            u_prov[i] = check;//big enough to see big errors if not handled correctly
        }
    }

    // if nn is not existent is being added in the boundary conditions part
    if (u_prov[0] != check && u_prov[1] != check){
        u[4] = 0.5 * (u_prov[0] + u_prov[1]);//v_n
    }
    else{
        u[4] = 0;
    }
    if (u_prov[2] != check && u_prov[3] != check){
        u[5] = 0.5 * (u_prov[2] + u_prov[3]);//v_s
    }
    else{
        u[5] = 0;
    }

    for (int i = 0; i < 4; ++i){
        if (u_prov[i + 4] != check){
            u[i] = 0.5 * (u_prov[i + 4] + u_P);
        }
        else{
            u[i] = 0;//no first neighbour in this direction, so no contribution
        }
    }

    double conv = - rho * ((u[4] * u[1] - u[5] * u[3]) * Sy + (u[0] * u[0] - u[2] * u[2]) * Sx);

    delete[] u;
    delete[] u_prov;

    return conv;
}
double conv_vv(double *u_n, double *v_n, double *pv_v, double rho, double Sx, double Sy, double v_P){//use CDS convective scheme
    double *v = new double[6]();// u_e, u_n, u_w, u_s, v_n, v_s
    double *v_prov = new double[8]();
    double check = 1e15;
    int pv_int;

    for (int i = 0; i < 4; ++i){//find velocities from u
        pv_int = (int)pv_v[i];
        if (pv_int != -1){
            v_prov[i] = u_n[pv_int];
        }
        else{
            v_prov[i] = check;//big enough to see big errors if not handled correctly
        }
    }
    for (int i = 4; i < 8; ++i){//find velocities from v
        pv_int = (int)pv_v[i];
        if (pv_int != -1){
            v_prov[i] = v_n[pv_int];
        }
        else{
            v_prov[i] = check;//big enough to see big errors if not handled correctly
        }
    }

    // if nn is not existent is being added in the boundary conditions part
    if (v_prov[0] != check && v_prov[3] != check){
        v[4] = 0.5 * (v_prov[0] + v_prov[3]);//u_e
    }
    else{
        v[4] = 0;
    }
    if (v_prov[1] != check && v_prov[2] != check){
        v[5] = 0.5 * (v_prov[1] + v_prov[2]);//u_w
    }
    else{
        v[5] = 0;
    }
    

    for (int i = 0; i < 4; ++i){
        if (v_prov[i + 4] != check){
            v[i] = 0.5 * (v_prov[i + 4] + v_P);
        }
        else{
            v[i] = 0;//no first neighbour in this direction, so no contribution
        }
    }

    double conv = - rho * ((v[1] * v[1] - v[3] * v[3]) * Sy + (v[4] * v[0] - v[5] * v[2]) * Sx);

    delete[] v;
    delete[] v_prov;

    return conv;
}

void boundary_U(int i, int nx, long int nn_u, double *u){
    if (i % (nx + 1) == 0){// left 
        u[i] = 0;
    }
    else if ((i + 1) % (nx + 1) == 0){//right
        u[i] = 0;
    }
    if(i + nx + 1 >= nn_u){//top
        u[i] = U_ref;
    }
    else if (i - nx - 1 < 0){//bottom
        u[i] = 0;
    }
}
void boundary_V(int i, int nx, long int nn_v, double *v){
    if ((i - nx < 0) || (i + nx) >= nn_v || (i % nx == 0) || (i + 1) % nx == 0){// bottom top left right
        v[i] = 0;
    } 
}
// double v_face_U(int j, int bnd_node, double dx, double dy, double val, double *pvi, double *v_n, int nn_u, int nx, double u_P){
//     double face_vel, S, v_ns;

//     if (j == 0){// right side
//         face_vel = u_P / 2; 
//         S = - dy;
//         return face_vel * face_vel * S;
//     }
//     else if (j == 1){// top side
//         if (bnd_node == nn_u - 1){//top right node
//             v_ns = 0; //mean velocity on the top right, find mean with 0 boundary cond
//         }
//         else if (bnd_node == nn_u - nx){//top left node
//             v_ns = 0; //mean velocity on the top left, find mean with 0 boundary cond
//         }
//         else{
//             v_ns = 0.5 * (v_n[(int)pvi[0]] + v_n[(int)pvi[1]]); //mean velocity on the top center
//         }
//         S = - dx;

//         face_vel = val;
//         return S * v_ns * face_vel;
//     }
//     else if (j == 2){// left side
//         face_vel = u_P / 2; 
//         S = dy;
//         return face_vel * face_vel * S;
//     }
//     else{// bottom side  if (j == 3)
//         if (bnd_node == nx - 1){//bottom right node
//             v_ns = 0; //mean velocity on the top right, find mean with 0 boundary cond
//         }
//         else if (bnd_node == 0){//bottom left node
//             v_ns = 0; //mean velocity on the top left, find mean with 0 boundary cond
//         }
//         else{
//             v_ns = 0.5 * (v_n[(int)pvi[2]] + v_n[(int)pvi[3]]); //mean velocity on the top center
//         }
//         S = dx;

//         face_vel = u_P / 2;
//         return S * v_ns * face_vel;
//     }
// }
// double v_face_V(int j, int bnd_node, double dx, double dy, double val, double *pvi, double *u_n, int nn_v, int nx, double v_P){
//     double face_vel, S, u_ew;

//     if (j == 0){// right side
//         if (bnd_node == nn_v - 1){//top right node
//             u_ew = 0; //mean velocity on the top right
//         }
//         else if (bnd_node == nx - 1){//bottom right node
//             u_ew = 0; //mean velocity on the bottom right
//         }
//         else{
//             u_ew = 0.5 * (u_n[(int)pvi[0]] + u_n[(int)pvi[3]]); //mean velocity on the right center
//         }
//         S = - dy;
//         face_vel = v_P / 2;

//         return S * u_ew * face_vel;
//     }
//     else if (j == 1){// top side
//         face_vel = v_P / 2; 
//         S = - dx;

//         return S * face_vel * face_vel;
//     }
//     else if (j == 2){// left side

//         if (bnd_node == nn_v - nx){//top left node
//             u_ew = 0; //mean velocity on the top left
//         }
//         else if (bnd_node == 0){//bottom left node
//             u_ew = 0; //mean velocity on the bottom left
//         }
//         else{
//             u_ew = 0.5 * (u_n[(int)pvi[1]] + u_n[(int)pvi[2]]); //mean velocity on the left center
//         }
//         S = dy;
//         face_vel = v_P / 2;

//         return S * u_ew * face_vel;
//     }
//     else{// bottom side  if (j == 3)
//         face_vel = v_P / 2; 
//         S = dx;

//         return S * face_vel * face_vel;
//     }
// }
void R_uv(double rho, int nx, int ny, double *u_n, double *R_u, double **pv_u, double **info_u, double **bound_u, long int nn_u, int nb_u, 
                                      double *v_n, double *R_v, double **pv_v, double **info_v, double **bound_v, long int nn_v, int nb_v){

    //u_n is the result at iteration n for u, R_u needs to be zeroed
    //v_n is the result at iteration n for v, R_v needs to be zeroed

    // allocate space for the two vectors, diffusion and convection for both v and u
    double *diff_u = new double[nn_u](), *conv_u = new double[nn_u]();
    double *diff_v = new double[nn_v](), *conv_v = new double[nn_v]();

    //definition of variables that will be used later
    double *tempI_U, *pvi_U, u_P, *bnd_U;
    double *tempI_V, *pvi_V, v_P, *bnd_V;

    double dx, dy;
    int bnd_node, nn = max(nn_u, nn_v), nb = max(nb_u, nb_v);

    #pragma omp parallel for num_threads(n_threads)
    for (int i = 0; i < nn; ++i){
        if (i < nn_u){
            tempI_U = info_u[i];// information matrix of point i
            dx = tempI_U[5];
            dy = tempI_U[6];
            pvi_U = pv_u[i];// first neighbours of point i
            u_P = u_n[i];// result of time n at point i

            // diffusion part
            R_u[i] = diff_uv(u_n, u_P, tempI_U, pvi_U);
        
            // convection part
            R_u[i] += conv_uu(u_n, v_n, pvi_U, rho, dy, dx, u_P);
        }
        if (i < nn_v){
            tempI_V = info_v[i];// information matrix of point i
            dx = tempI_V[5];
            dy = tempI_V[6];
            pvi_V = pv_v[i];// first neighbours of point i
            v_P = v_n[i];// result of time n at point i
    
            // diffusion part
            R_v[i] = diff_uv(v_n, v_P, tempI_V, pvi_V);
        
            // convection part
            R_v[i] += conv_vv(u_n, v_n, pvi_V, rho, dy, dx, v_P);
        }
    }

    // double face_vel, S;

    // //boundary conditions
    // #pragma omp parallel for num_threads(n_threads)
    // for (int k = 0; k < nb; ++k){
    //     if (k < nb_u){
    //         bnd_U = bound_u[k];
    //         bnd_node = (int)bnd_U[0];
    //         pvi_U = pv_u[bnd_node];//first neighbours
    //         tempI_U = info_u[bnd_node];// information matrix of point bnd_node
    //         dx = tempI_U[5];
    //         dy = tempI_U[6];
    //         u_P = u_n[bnd_node];// result of time n at point bnd_node
            
    //         int j = (int)bnd_U[3];
            
    //         // keep only dirichlet boundary conditions as only those are present in the problem to be solved
    //         if (bnd_U[2] == 2){//pure dirichlet non time dependent
    //             diff_u[bnd_node] += tempI_U[j + 1] * (bnd_U[1] - u_P);
    //             conv_u[bnd_node] += rho * v_face_U(j,bnd_node,dx,dy,bnd_U[1],pvi_U,v_n,nn_u,nx+1,u_P);
    //         }
    //     }
    //     if (k < nb_v){
    //         bnd_V = bound_v[k];
    //         bnd_node = (int)bnd_V[0];
    //         pvi_V = pv_v[bnd_node];//first neighbours
    //         tempI_V = info_v[bnd_node];// information matrix of point bnd_node
    //         dx = tempI_V[5];
    //         dy = tempI_V[6];
    //         v_P = v_n[bnd_node];// result of time n at point bnd_node

    //         int j = (int)bnd_V[3];
            
    //         // keep only dirichlet boundary conditions as only those are present in the problem to be solved
    //         if (bnd_V[2] == 2){//pure dirichlet non time dependent
    //             diff_v[bnd_node] += tempI_V[j + 1] * (bnd_V[1] - v_P);
    //             conv_v[bnd_node] += rho * v_face_V(j,bnd_node,dx,dy,bnd_V[1],pvi_V,u_n,nn_v,nx,v_P);
    //         }
    //     }
    // }

    // vector_up(diff_u, conv_u, R_u, nn_u, 0); //conv+diff for u
    // vector_up(diff_v, conv_v, R_v, nn_v, 0); //conv+diff for v

    // Clean up dynamically allocated memory
    delete[] diff_u;
    delete[] conv_u;
    delete[] diff_v;
    delete[] conv_v;
}
// void pre_bound(int i, int j, int nx, int ny, double *vel, double **info_uv, double *dd, double *uv, int direction){
//     // vel must be passed as a vector of length 8, containing in the values for the convection and the diffusion (ENWS,enws)
//     // uv is the vector of the velocities, in u or in v depending which is passed
//     // pvi is the first nearest neighbours of node ij

//     long int node = i + nx * j, nodeE = node + 1, nodeW = node - 1, nodeN = node + nx, nodeS = node - nx;
//     double uv_P = uv[node];

//     if (i < nx - 1){// not right
//         vel[0] = uv[nodeE];
//         dd[0] = (info_uv[node][5] + info_uv[nodeE][5]) / 2;
//     }
//     else{
//         vel[0] = 0;
//         dd[0] = info_uv[node][5];
//     }
//     if (j < ny - 1){//not top
//         vel[1] = uv[nodeN];
//         dd[1] = (info_uv[node][6] + info_uv[nodeN][6]) / 2;
//     }
//     else{
//         if (direction == 0){//u
//             vel[1] = 2 * U_ref - uv_P;
//             dd[1] = info_uv[node][6];
//         }
//         else{//v
//             vel[1] = 0;
//             dd[1] = info_uv[node][6];
//         }
//     }
//     if (i > 0){//not left
//         vel[2] = uv[nodeW];
//         dd[2] = (info_uv[node][5] + info_uv[nodeW][5]) / 2;
//     }
//     else{
//         vel[2] = 0;
//         dd[2] = info_uv[node][5];
//     }
//     if (j > 0){// not bottom
//         vel[3] = uv[nodeS];
//         dd[3] = (info_uv[node][6] + info_uv[nodeS][6]) / 2;
//     }
//     else{
//         vel[3] = 0;
//         dd[3] = info_uv[node][6];
//     }

//     for (int i = 4; i < 8; ++i){
//         vel[i] = 0.5 * (vel[i - 4] + uv_P);
//     }
// }
// void pre_bound_vel1(int i, int j, int nx, int ny, double *vel, double *uv, int direction){
//     long int node = i + nx * j, nodeE = node + 1, nodeW = node - 1, nodeN = node + nx, nodeS = node - nx;
    
//     if (direction == 0){//u
//         if (i == 0 || i == nx - 1){
//             vel[0] = 0;
//             vel[1] = 0;
//         }
//         else{
//             vel[0] = 0.5 * (uv[nodeN - 1] + uv[nodeN]);
//             vel[1] = 0.5 * (uv[nodeW] + uv[node]);    
//         }
//     }
//     else{
//         if (j == 0 || j == ny - 1){
//             vel[0] = 0;
//             vel[1] = 0;
//         }
//         else{
//             vel[0] = 0.5 * (uv[node] + uv[nodeS]);
//             vel[1] = 0.5 * (uv[nodeE] + uv[nodeS + 1]);    
//         }
//     }
// }
// double diffus(double dx, double dy, double *vel, double uv_P, double *dist){
//     double diff = 0;
//     double S;

//     for (int k = 0; k < 4; ++k){
//         if (k == 0 || k == 2){
//             S = dy;
//         }
//         else{
//             S = dx;
//         }
//         diff += (vel[k] - uv_P) * S / dist[k];//put correct distance
//     }

//     return diff;
// }
// void R_uv(double rho, double mu, int nx, int ny, double *u_n, double *R_u, double **info_u, long int nn_u, 
//                                                  double *v_n, double *R_v, double **info_v, long int nn_v){

//     double *diff_u = new double[nn_u](), *conv_u = new double[nn_u]();
//     double *diff_v = new double[nn_v](), *conv_v = new double[nn_v]();
//     double *vel = new double[8]();
//     double *vel1 = new double[2]();
//     double *dist = new double[4]();

//     //definition of variables that will be used later
//     double *tempI_U, *pvi_U, u_P;
//     double *tempI_V, *pvi_V, v_P;

//     double dx, dy;
//     long int node;

//     for (int j = 0; j < ny + 1; ++j){
//         for (int i = 0; i < nx + 1; ++i){
//             if (j < ny){//u
//                 node = i + (nx + 1) * j;
//                 tempI_U = info_u[node];// information matrix of point i
//                 dx = tempI_U[5];
//                 dy = tempI_U[6];

//                 pre_bound(i, j, nx+1, ny, vel, info_u, dist, u_n, 0);
//                 pre_bound_vel1(i, j, nx, ny, vel1, v_n, 0);
//                 conv_u[node] = - rho * (vel[4]*vel[4]*dy + vel1[0]*vel[5]*dx - vel[6]*vel[6]*dy - vel1[1]*vel[7]*dx);
                
//                 diff_u[node] = mu * diffus(dx,dy,vel,u_n[node],dist);
//             }
//             if (i < nx){//v
//                 node = i + nx * j;
//                 tempI_V = info_v[node];// information matrix of point i
//                 dx = tempI_V[5];
//                 dy = tempI_V[6];
//                 pre_bound(i, j, nx, ny+1, vel, info_v, dist, v_n, 1);
//                 pre_bound_vel1(i, j, nx+1, ny, vel1, u_n, 1);
//                 conv_v[node] = - rho * (vel1[1]*vel[4]*dy + vel[5]*vel[5]*dx - vel1[0]*vel[6]*dy - vel[7]*vel[7]*dx);
                
//                 diff_v[node] = mu * diffus(dx,dy,vel,v_n[node],dist);
//             }
//         }
//     }

//     vector_up(diff_u, conv_u, R_u, nn_u, 0); //conv+diff for u
//     vector_up(diff_v, conv_v, R_v, nn_v, 0); //conv+diff for v

//     // Clean up dynamically allocated memory
//     delete[] diff_u;
//     delete[] conv_u;
//     delete[] diff_v;
//     delete[] conv_v;
//     delete[] vel;
//     delete[] vel1;
//     delete[] dist;
// }
void compute_R_uv_eff(double **R_u, double *R_uEff, double **R_v, double *R_vEff, long int nn_u, long int nn_v){
    long int nn = max(nn_u,nn_v);
    double *RU0 = R_u[0], *RU1 = R_u[1];
    double *RV0 = R_v[0], *RV1 = R_v[1];

    #pragma omp parallel for num_threads(n_threads)
    for (int i = 0; i < nn; ++i){
        if (i < nn_u){
            R_uEff[i] = 1.5 * RU1[i] - 0.5 * RU0[i];
        }
        if (i < nn_v){
            R_vEff[i] = 1.5 * RV1[i] - 0.5 * RV0[i];
        }
    }
}
void compute_trial_vel(double dt, double rho, int nx, double **info_u, long int nn_u, double *u_trial, double *u_n, double *R_u,
                                                      double **info_v, long int nn_v, double *v_trial, double *v_n, double *R_v){
    //don't need trial velocity at previous time step so i zero to have a blank vector
    // zero(u_trial,nn_u);
    // zero(v_trial,nn_v);

    int nn = max(nn_u, nn_v);
    double K, *tempIU, *tempIV;

    #pragma omp parallel for num_threads(n_threads)
    for (int i = 0; i < nn; ++i){
        if (i < nn_u){
            tempIU = info_u[i];
            K = dt / (rho * tempIU[0]);// constant to multiply R_uv
            u_trial[i] = u_n[i] + K * R_u[i];

            boundary_U(i,nx,nn_u,u_trial);
        }
        if (i < nn_v){
            tempIV = info_v[i];
            K = dt / (rho * tempIV[0]);// constant to multiply R_uv
            v_trial[i] = v_n[i] + K * R_v[i];

            boundary_V(i,nx,nn_v,v_trial);
        }
    }
}

// void new_vel(double dt_rho, double *P_n1, long int nn_p, int nx, double **info_u, double *u_trial, double *u_n1, long int nn_u, 
//                                                                  double **info_v, double *v_trial, double *v_n1, long int nn_v){//dt_rho = dt/rho

//     double dx, dy, *tempIU, *tempIV;
//     int contP_u = 0, contP_v = nx;//start counting the pressure point from -1 as I need the first point to the left to get the pressure in point 0 for u
//     long int nn = max(nn_u, nn_v);

//     //cycle over all the nodes
//     for (int i = 0; i < nn; ++i){
//         if (i < nn_u){
//             tempIU = info_u[i];
//             dx = tempIU[5];
//             dy = tempIU[6];
//             u_n1[i] = u_trial[i] - (dt_rho / dx) * (P_n1[contP_u + 1] - P_n1[contP_u]);

//             // boundary_U(i,nx,nn_u,u_n1);
//             if ((i - nx) % (nx + 1) == 0 || i == 0 || contP_u == nn_p - 2){//right, the other two conditions to not make p go out of bounds
//                 contP_u--;//to make correct the next pressure point
//             }

//             contP_u++;//next pressure point
//         }
//         if (i < nn_v){
//             tempIV = info_v[i];
//             dx = tempIV[5];
//             dy = tempIV[6];
//             v_n1[i] = v_trial[i] - (dt_rho / dy) * (P_n1[contP_v] - P_n1[contP_v - nx]);

//             // boundary_V(i,nx,nn_v,v_n1);
//             if (i < nx + 1){
//                 contP_v--;
//             }

//             contP_v++;//next pressure point
//         }
//     }
// }

void new_vel(double dt_rho, double *P_n1, int nx, int ny, double **info_u, double *u_trial, double *u_n1, 
                                                          double **info_v, double *v_trial, double *v_n1){//dt_rho = dt/rho

    double dx, dy, *tempIU, *tempIV;
    int contP_u = -1, contP_v = 0;//start counting the pressure point from -1 as I need the first point to the left to get the pressure in point 0 for u
    long int nn_u = ny * (nx + 1), nn_v = nx * (ny + 1), node, P_next, P_prev;


    // copy_vector(u_n1, u_trial, nn_u);
    // copy_vector(v_n1, v_trial, nn_v);

    for (int j = 0; j < ny + 1; ++j){
        for (int i = 0; i < nx + 1; ++i){
            if(j < ny){
                node = i + (nx + 1) * j;
                u_n1[node] = u_trial[node];

                if (i > 0 && i < nx){//u
                    P_next = i + nx * j;
                    P_prev = i - 1 + nx * j;
                    tempIU = info_u[node];
                    dx = tempIU[5];
                    dy = tempIU[6];
        
                    u_n1[node] -= (dt_rho / dx) * (P_n1[P_next] - P_n1[P_prev]);
                    // boundary_U(i,nx,nn_u,u_n1);
                }
            }
            if(i < nx){
                node = i + nx * j;
                v_n1[node] = v_trial[node];

                if (j > 0 && j < ny){//v
                    P_next = i + nx * j;
                    P_prev = i + nx * (j - 1);
                    tempIV = info_v[node];
                    dx = tempIV[5];
                    dy = tempIV[6];
        
                    v_n1[node] -= (dt_rho / dy) * (P_n1[P_next] - P_n1[P_prev]);
                    // boundary_V(i,nx,nn_v,v_n1);
                }
            }
        }
    }
}
double *vel_mod(double *u_fin, double *v_fin, long int nn_p){
    int contx = 0, nx = sqrt(nn_p);
    double *v_mod = new double[nn_p]();

    // #pragma omp parallel for num_threads(n_threads)
    for (int i = 0; i < nn_p; ++i){
        if ((i - nx) % (nx + 1) == 0){//right
            contx++;
            v_mod[i] = sqrt(u_fin[contx]*u_fin[contx] + v_fin[i]*v_fin[i]);
        }
        else{
            v_mod[i] = sqrt(u_fin[contx]*u_fin[contx] + v_fin[i]*v_fin[i]);
        }
        contx++;
    }
    return v_mod;
}














#endif
