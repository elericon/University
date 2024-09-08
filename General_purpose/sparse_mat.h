#ifndef SPARSE_MAT_H
#define SPARSE_MAT_H



class sparse{
	private:
		int n_row,n_term,n_col;//n_col not really needed but useful info for checks in case of need
		double *v_coef;//value 
		int *Ja;//column indices
		int *iat;//pointer to first nonzero entry of each row
	public: 
	
		//constructor and destructor
		sparse();
		sparse(int nn, int nt);
		sparse(int nr, int nc, int nt);
		~sparse();

		//getters 
		int get_nrow(){return n_row;}
		int get_ncol(){return n_col;}
		int get_nterm(){return n_term;}
		double *get_coef(){return v_coef;}
		int *get_Ja(){return Ja;}
		int *get_iat(){return iat;}

		//setter
		void set_nrow(int n){n_row = n;}
		void set_ncol(int n){n_col = n;}
		void set_nterm(int n){n_term = n;}
		void set_coef(double *v){v_coef = v;}
		void set_Ja(int *v){Ja = v;}
		void set_iat(int *v){iat = v;}

		void alloc(int nr, int nt){
			n_row = nr;
			n_term = nt;
			n_col = nr;
			v_coef = new double[nt]();
			Ja = new int[nt]();
			iat = new int[nr+1]();
		}
		void alloc_one(int option, int n){
			if(option == 0){
				n_row = n;
			}
			else if(option == 1){
				n_term = n;
			}
			else if(option == 2){
				v_coef = new double[n]();
			}
			else if(option == 3){
				Ja = new int[n]();
			}
			else if(option == -1){
				n_col = n;
			}
			else{
				iat = new int[n]();
			}
		}
		
		void copy_Ja(int *v){
			#pragma omp parallel for num_threads(n_threads)
			for (int i = 0; i < n_term; ++i){
				Ja[i] = v[i];
			}
		}
		void copy_iat(int *v){
			#pragma omp parallel for num_threads(n_threads)
			for (int i = 0; i < n_row+1; ++i){
				iat[i] = v[i];
			}
		}
		void copy_coef(double *v){
			#pragma omp parallel for num_threads(n_threads)
			for (int i = 0; i < n_term; ++i){
				v_coef[i] = v[i];
			}
		}
		void copy(sparse &original){//the original matrix needs to have the same space allocated for the one calling this function
			n_col = original.get_ncol();
			double *v_O = original.get_coef();
			int *iat_O = original.get_iat(), *ja_O = original.get_Ja();

			#pragma omp parallel for num_threads(n_threads)
			for (int i = 0; i < n_term; ++i){
				Ja[i] = ja_O[i];
				v_coef[i] = v_O[i];
			}
			#pragma omp parallel for num_threads(n_threads)
			for (int i = 0; i < n_row+1; ++i){
				iat[i] = iat_O[i];
			}
		}
		void print(){
			std::cout<<"Vector of coefficients:\n";
			print_vett(v_coef,n_term);
			std::cout<<"Vector Ja:\n";
			print_vett(Ja,n_term);
			std::cout<<"Vector Iat:\n";
			print_vett(iat,n_row+1);
		}
		void print(std::string filename){
			std::ofstream print;
			print.open(filename);
			print << n_row << "\n";
			print << n_term << "\n";
			for (int i = 0; i < n_term; ++i){
				print << v_coef[i] << "\t";
			}
			print << "\n";
			for (int i = 0; i < n_term; ++i){
				print << Ja[i] << "\t";
			}
			print << "\n";
			for (int i = 0; i < n_row+1; ++i){
				print << iat[i] << "\t";
			}
			print << "\n";
			print.close();
		}
		void import(std::string filename){
			std::ifstream print;
			print.open(filename);
			print >> n_row;
			print >> n_term;
			alloc(n_row,n_term);
			for (int i = 0; i < n_term; ++i){
				print >> v_coef[i];
			}
			for (int i = 0; i < n_term; ++i){
				print >> Ja[i];
			}
			for (int i = 0; i < n_row+1; ++i){
				print >> iat[i];
			}
			print.close();
		}
		void sumv(double *v){//assume i'm adding a matrix with the same sparsity pattern 
			#pragma omp parallel for num_threads(n_threads)
			for (int i = 0; i < n_term; ++i){
				v_coef[i] += v[i];
			}
		}
		void subv(double *v){//assume i'm subtracting a matrix with the same sparsity pattern 
			#pragma omp parallel for num_threads(n_threads)
			for (int i = 0; i < n_term; ++i){
				v_coef[i] -= v[i];
			}
		}
		void multiply(double k){
			#pragma omp parallel for num_threads(n_threads)
			for (int i = 0; i < n_term; ++i){
				v_coef[i] *= k;
			}
		}
		void divide(double k){
			#pragma omp parallel for num_threads(n_threads)
			for (int i = 0; i < n_term; ++i){
				v_coef[i] /= k;
			}
		}
		sparse diag(){
			// allocate space for the diagonal matrix 
			double *v_coefj = new double[n_row];
			int *Jaj = new int[n_row],*iatj = new int[n_row+1];
		
			int iaend = iat[0], iastart;

			#pragma omp parallel for private(iastart,iaend) num_threads(n_threads)
			for (int i = 0; i < n_row; ++i){//cycle over iat
				iastart = iat[0];
				iaend = iat[i+1];
				for (int j = iastart; j < iaend; ++j){//cycle on part of iat
					if(Ja[j] == i){
						v_coefj[i] = v_coef[j];
					}
					Jaj[i] = i;
					iatj[i] = i;
				}
			}
			iatj[n_row] = n_row;
		
			sparse jac(n_row, n_row);
			jac.copy_coef(v_coefj);
			jac.copy_iat(iatj);
			jac.copy_Ja(Jaj);

			return jac;
		}
		sparse diag_jac(){
			// allocate space for the jacobi vectors, remember diagonal matrix (computed as the preconditioner M)
			double *v_coefj = new double[n_row];
			int *Jaj = new int[n_row],*iatj = new int[n_row+1];
		
			int iaend = iat[0], iastart;

			#pragma omp parallel for private(iastart,iaend) num_threads(n_threads)
			for (int i = 0; i < n_row; ++i){//cycle over iat
				iastart = iat[i];
				iaend = iat[i+1];
				for (int j = iastart; j < iaend; ++j){//cycle on part of iat
					if(Ja[j] == i){
						if (v_coef[j] != 0){
							v_coefj[i] = 1/v_coef[j];
						}
						else{
							v_coefj[i] = 1;
						}
					}
					Jaj[i] = i;
					iatj[i] = i;
				}
			}
			iatj[n_row] = n_row;
		
			sparse jac(n_row, n_row);
			jac.copy_coef(v_coefj);
			jac.copy_iat(iatj);
			jac.copy_Ja(Jaj);

			return jac;
		}
		// creates a sparse matrix which is equivalent to [A, B] in MATLAB, need space already allocated
		void sparse_assemble_mat_row(sparse &A, sparse &B){
			int n_colA = A.get_ncol();
			n_col = n_colA + B.get_ncol();
			// get handles to speed up
			double *vA = A.get_coef();
			int *jaA = A.get_Ja(), *iatA = A.get_iat();
		
			// get handles to speed up
			double *vB = B.get_coef();
			int *jaB = B.get_Ja(), *iatB = B.get_iat();
		
			// iat just need to sum the two separate iat to become the new one
			vector_up(iatA, iatB, iat, n_row + 1, 0);
		
			// define start and end of the for cycles
			int iastartA = iatA[0], iastartB = iatB[0], iastartM = iat[0];
		
			// cycle across the iatM vector
			for (int cont = 1; cont < n_row + 1; ++cont){
				// choose ending points of iatABM
				int iaendA = iatA[cont], iaendB = iatB[cont], iaendM = iat[cont], endA = iaendA - iastartA;
		
				// add the jaA, vA to jaM and vM
				for (int i = 0; i < endA; ++i){
					Ja[i + iastartM] = jaA[i + iastartA];
					v_coef[i + iastartM] = vA[i + iastartA];
				}
		
				// add the jaB, vB to jaM and vM
				for (int i = 0; i < iaendB - iastartB; ++i){
					Ja[i + iastartM + endA] = jaB[i + iastartB] + n_colA;
					v_coef[i + iastartM + endA] = vB[i + iastartB];
				}
		
				// choose next starting points of iatABM
				iastartA = iaendA;
				iastartB = iaendB;
				iastartM = iaendM;
			}
		}
		// creates a sparse matrix which is equivalent to [A; B] in MATLAB, need space already allocated
		void sparse_assemble_mat_col(sparse &A, sparse &B){
			int ntA = A.get_nterm(), ntB = B.get_nterm();
			int nrA = A.get_nrow(), nrB = B.get_nrow();

			// define handles for speed
			double *vA = A.get_coef(), *vB = B.get_coef();
			int *iatA = A.get_iat(), *iatB = B.get_iat(), *jaA = A.get_Ja(), *jaB = B.get_Ja();

			//fix the coef and ja for A
			for (int i = 0; i < ntA; ++i){
				v_coef[i] = vA[i];
				Ja[i] = jaA[i];
			}
			//fix the coef and ja for B
			for (int i = 0; i < ntB; ++i){
				v_coef[i + ntA] = vB[i];
				Ja[i + ntA] = jaB[i];
			}

			//fix iat for A
			for (int i = 0; i < nrA; ++i){
				iat[i] = iatA[i];
			}
			//fix iat for B
			int endA = iatA[nrA];
			for (int i = 0; i < nrB + 1; ++i){
				iat[i + nrA] = iatB[i] + endA;
			}
		}
		// creates a sparse matrix which is equivalent to [A, 0; 0, B] in MATLAB, need space already allocated
		void sparse_assemble_blockdiag(sparse &A, sparse &B){
			int ntA = A.get_nterm(), ntB = B.get_nterm();
			int nrA = A.get_nrow(), nrB = B.get_nrow();

			// define handles for speed
			double *vA = A.get_coef(), *vB = B.get_coef();
			int *iatA = A.get_iat(), *iatB = B.get_iat(), *jaA = A.get_Ja(), *jaB = B.get_Ja();

			//fix the coef and ja for A
			for (int i = 0; i < ntA; ++i){
				v_coef[i] = vA[i];
				Ja[i] = jaA[i];
			}
			//fix the coef and ja for B
			for (int i = 0; i < ntB; ++i){
				v_coef[i + ntA] = vB[i];
				Ja[i + ntA] = jaB[i] + nrA;// traslate the values by a nrB x nrA matrix of zeros
			}

			//fix iat for A
			for (int i = 0; i < nrA; ++i){
				iat[i] = iatA[i];
			}
			//fix iat for B
			int endA = iatA[nrA];
			for (int i = 0; i < nrB + 1; ++i){
				iat[i + nrA] = iatB[i] + endA;
			}
		}
		void transpose(sparse &original) {//number of columns is the number of rows of the matrix using this function. space must already be allocated.
		    int nr_orig = original.get_nrow();
		
		    // Temporary storage for the number of entries in each column of the original matrix
		    int *col_counts = new int[n_row]();

		    //create handles of the original matrix
		    int *ja_tbt = original.get_Ja();
		    int *iat_tbt = original.get_iat();
		    double *v_tbt = original.get_coef();
		
		    // Count the number of entries in each column (these will become the rows of the transposed matrix)
		    for (int i = 0; i < n_term; ++i) {
		        col_counts[ja_tbt[i]]++;
		    }
		
		    // Set the iat array for the transposed matrix
		    iat[0] = 0;
		    for (int i = 0; i < n_row; ++i) {
		        iat[i + 1] = iat[i] + col_counts[i];
		    }
		
		    // Temporary storage for the current position in each column
		    int *current_pos = new int[n_row]();
		    for (int i = 0; i < n_row; ++i) {
		        current_pos[i] = iat[i];
		    }
		
		    // Fill the transposed matrix
		    for (int i = 0; i < nr_orig; ++i) {
		        for (int j = iat_tbt[i]; j < iat_tbt[i + 1]; ++j) {
		            int col = ja_tbt[j];
		            int pos = current_pos[col];
		            Ja[pos] = i;
		            v_coef[pos] = v_tbt[j];
		            current_pos[col]++;// move onto iat new to get the correct space in which to fill
		        }
		    }
		
		    // Clean up temporary arrays
		    delete[] col_counts;
		    delete[] current_pos;
		}
		void zeros(){
			#pragma omp parallel for num_threads(n_threads)
			for (int i = 0; i < n_term; ++i){
				v_coef[i] = 0;
			}
		}
		void zero_all(){
			#pragma omp parallel for num_threads(n_threads)
			for (int i = 0; i < n_term; ++i){ 
				v_coef[i] = 0;
				Ja[i] = 0;
			}
			#pragma omp parallel for num_threads(n_threads)
			for(int i = 0; i < n_row; ++i){
				iat[i] = 0;
			}
		}

		void id(){//assume the sparsity pattern permits nonzero diagonal
			for (int i = 0; i < n_row + 1; ++i){
				for (int j = iat[i]; j < iat[i+1]; ++j){
					if (Ja[j] == i){
						v_coef[j] = 1;
						break;
					}
				}
			}
		}
};





sparse :: sparse(int nn, int nt){//assume square
	n_row = nn;
	n_term = nt;
	n_col = nn;
	v_coef = ::new double[nt]();
	Ja = ::new int[nt]();
	iat = ::new int[nn+1]();
};
sparse :: sparse(int nr, int nc, int nt){
	n_row = nr;
	n_term = nt;
	n_col = nc;
	v_coef = ::new double[nt]();
	Ja = ::new int[nt]();
	iat = ::new int[nr+1]();
}
sparse :: sparse(){
}
sparse :: ~sparse(){
	delete[] v_coef;
	delete[] Ja;
	delete[] iat;
};

// Needs mat to have space allocated already
void render_sparse(double **M, sparse *mat_ptr){//suppose the matrix is a square matrix and compute it's sparse counterpart
    sparse &mat = *mat_ptr;

    int row = mat.get_nrow();
    int n_term = mat.get_nterm();

    // Define the handles
    double *v_coef = mat.get_coef();
    int *Ja = mat.get_Ja();
    int *iat = mat.get_iat();
    double *row_i, temp;

    // Assign the nonzero values to the sparse matrix
    int cont = 0; // For the coefficient
    
   
    // #pragmma omp parallel for schedule(dynamic, row/1000+1)
    for (int i = 0; i < row; ++i) {
        iat[i] = cont; // Update iat for the current row
        row_i = M[i];

        for (int j = 0; j < row; ++j) {
        	temp = row_i[j];//only access once

            if (temp != 0) {
                v_coef[cont] = temp;
                Ja[cont] = j;
               
                // #pragmma omp atomic
                cont++;
            }
        }
    }
    iat[row] = n_term; // Update the last entry in iat
}
void render_sparse(double **M, sparse *mat_ptr, int nc) {
    sparse &mat = *mat_ptr;

    int row = mat.get_nrow();
    int n_term = mat.get_nterm();

    // Define the handles
    double *v_coef = mat.get_coef();
    int *Ja = mat.get_Ja(), *iat = mat.get_iat();
    double *row_i, temp;

    // Assign the nonzero values to the sparse matrix
    int cont = 0; // For the coefficient

    for (int i = 0; i < row; ++i) {
        iat[i] = cont; // Update iat for the current row
        row_i = M[i];

        for (int j = 0; j < nc; ++j) {
            temp = row_i[j]; // Only access once

            if (temp != 0) {
                v_coef[cont] = temp;
                Ja[cont] = j;
                cont++;
            }
        }
    }
    iat[row] = n_term; // Update the last entry in iat
}

void render_full(double **M, sparse &mat){//the full matrix needs to have space already allocated
	//define handles
	int row = mat.get_nrow();
	int nterm = mat.get_nterm();
	double *v_coef = mat.get_coef();
	int *Ja = mat.get_Ja();
	int *iat = mat.get_iat();

	//define the counters
	int cont1 = 0;

	//complete matrix M
	int jend = iat[0], jstart;
	// #pragmma omp parallel for schedule(dynamic,row/1000+1) private(jend,jstart)
	for (int i = 0; i < row; ++i){
		for (int j = 0; j < row; ++j){//sets all the row to zeroes
			M[i][j] = 0;
		}
		jstart = iat[i];
		jend = iat[i+1];
		
		for (int j = jstart; j < jend; ++j){//modify only the entries
			M[i][Ja[j]] = v_coef[cont1];
			// #pragmma omp atomic
			cont1++;
		}
	}
}
void sparse_vm_prod(sparse &mat, double *x, double *b){// computes sparse matrix vector product 
	// handles
	int n_row = mat.get_nrow(), n_term = mat.get_nterm();
	double *v_coef;
	int *Ja,*iat;
	v_coef = mat.get_coef();
	Ja = mat.get_Ja();
	iat = mat.get_iat();


	//scalar product
	int iend = iat[0];
	#pragma omp parallel for private(iend) num_threads(n_threads)
	for (int i = 0; i < n_row; ++i){
		b[i] = 0;
		int istart = iat[i];
		iend = iat[i+1];
		for (int j = istart; j < iend; ++j){
			b[i] += v_coef[j]*x[Ja[j]];
		}
	}
}
void copy_spMat(sparse &result, sparse &to_copy){//suppose already allocated the space
	result.copy_Ja(to_copy.get_Ja());
	result.copy_iat(to_copy.get_iat());
	result.copy_coef(to_copy.get_coef());
}

// void sparse_sum(sparse &A, sparse &B) {
//   // Check if matrices have compatible dimensions
//   if (A.get_nrow() != B.get_nrow() || A.get_ncol() != B.get_ncol()) {
//     throw std::invalid_argument("Matrices must have the same dimensions for summation");
//   }

//   // Allocate space for the resulting sparse matrix (assume worst-case scenario: all elements non-zero)
//   int n_row = A.get_nrow();
//   int n_col = A.get_ncol();
//   int max_nt = A.get_nterm() + B.get_nterm();
//   sparse C(n_row, n_col, max_nt);

//   // Initialize counters for non-zero elements in A, B, and C
//   int i_A = 0, i_B = 0, i_C = 0;

//   // Loop through each row of the matrices
//   for (int i = 0; i < n_row; ++i) {
//     // Loop through elements in row i of A
//     while (i_A < A.iat[i + 1] && A.Ja[i_A] < n_col) {
//       int col_A = A.Ja[i_A];

//       // Loop through elements in row i of B
//       while (i_B < B.iat[i + 1] && B.Ja[i_B] < n_col) {
//         int col_B = B.Ja[i_B];

//         // If both A and B have elements in the same column
//         if (col_A == col_B) {
//           C.v_coef[i_C] = A.v_coef[i_A] + B.v_coef[i_B];
//           C.Ja[i_C] = col_A;
//           i_A++;
//           i_B++;
//           i_C++;
//         } else if (col_A < col_B) {
//           // Copy element from A (no corresponding element in B)
//           C.v_coef[i_C] = A.v_coef[i_A];
//           C.Ja[i_C] = col_A;
//           i_A++;
//           i_C++;
//         } else {
//           // Copy element from B (no corresponding element in A)
//           C.v_coef[i_C] = B.v_coef[i_B];
//           C.Ja[i_C] = col_B;
//           i_B++;
//           i_C++;
//         }
//       }

//       // Check if there are remaining elements in A for this row
//       if (i_A < A.iat[i + 1]) {
//         // Copy remaining elements from A (no corresponding elements in B)
//         C.v_coef[i_C] = A.v_coef[i_A];
//         C.Ja[i_C] = A.Ja[i_A];
//         i_A++;
//         i_C++;
//       }
//       // Move to the next element in B
//       i_B = B.iat[i];
//     }
//   }

//   // Update the number of non-zero elements in the resulting matrix
//   C.set_nterm(i_C);

//   // No need to modify iat since it's already populated during the loop
// }

#endif
