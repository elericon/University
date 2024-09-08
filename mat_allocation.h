#ifndef MAT_ALLOCATION_H
#define MAT_ALLOCATION_H


double **mat_alloc(int nr, int nc){//returns the allocated matrix as pointer to pointer and already initialized to zero
	if(nr == 0){
		std::cout<<"Error, the matrix has row dimension zero, call: "<<numcall<<"\n";
		return nullptr;
	}
	if(nc == 0){
		std::cout<<"Error, the matrix has col dimension zero, call: "<<numcall<<"\n";
		return nullptr;
	}
	
	numcall++;
		
	double **M = new double*[nr]();
	double *buffer = new double[nr*nc]();

	for (int i = 0; i < nr; ++i){
		M[i] = &buffer[i*nc];
	}
	// zero(M,nr,nc);
	return M;
}
double **mat_alloc(long int nr, long int nc){//returns the allocated matrix as pointer to pointer and already initialized to zero
	if(nr == 0){
		std::cout<<"Error, the matrix has row dimension zero, call: "<<numcall<<"\n";
		return nullptr;
	}
	if(nc == 0){
		std::cout<<"Error, the matrix has col dimension zero, call: "<<numcall<<"\n";
		return nullptr;
	}
	
	numcall++;
		
	double **M = new double*[nr]();
	double *buffer = new double[nr*nc]();

	for (int i = 0; i < nr; ++i){
		M[i] = &buffer[i*nc];
	}
	// zero(M,nr,nc);
	return M;
}
double **mat_alloc(long int nr, int nc){//returns the allocated matrix as pointer to pointer and already initialized to zero
	if(nr == 0){
		std::cout<<"Error, the matrix has row dimension zero, call: "<<numcall<<"\n";
		return nullptr;
	}
	if(nc == 0){
		std::cout<<"Error, the matrix has col dimension zero, call: "<<numcall<<"\n";
		return nullptr;
	}
	
	numcall++;
		
	double **M = new double*[nr]();
	double *buffer = new double[nr*nc]();

	for (int i = 0; i < nr; ++i){
		M[i] = &buffer[i*nc];
	}
	// zero(M,nr,nc);
	return M;
}
double **mat_alloc(int nr, long int nc){//returns the allocated matrix as pointer to pointer and already initialized to zero
	if(nr == 0){
		std::cout<<"Error, the matrix has row dimension zero, call: "<<numcall<<"\n";
		return nullptr;
	}
	if(nc == 0){
		std::cout<<"Error, the matrix has col dimension zero, call: "<<numcall<<"\n";
		return nullptr;
	}
	
	numcall++;
		
	double **M = new double*[nr]();
	double *buffer = new double[nr*nc]();

	for (int i = 0; i < nr; ++i){
		M[i] = &buffer[i*nc];
	}
	// zero(M,nr,nc);
	return M;
}

void mat_dealloc(double **M){
	delete[] M[0];
	delete[] M;
}
double **assemble_mat(double **A, double **B, double **C, int nr, int nc){
	double **M = mat_alloc(nr + nc, nr + nc);
	double *temp,*temp1;
	for (int i = 0; i < nr; ++i){
		temp = M[i];
		temp1 = A[i];
		for (int j = 0; j < nr; ++j){
			temp[j] = temp1[j];
		}
		for (int j = 0; j < nc; ++j){
			temp[j + nr] = B[j][i];
		}
	}
	for (int i = 0; i < nc; ++i){
		temp = M[i + nr];
		temp1 = B[i];
		for (int j = 0; j < nr; ++j){
			temp[j] = temp1[j];
		}
		temp1 = C[i];
		for (int j = 0; j < nc; ++j){
			temp[j + nr] = temp1[j];
		}
	}
	return M;
}
double **import_mat(std::string filename, int &nr, int &nc){
	double line;
	std::ifstream file;
	file.open(filename);
	file >> nr >> nc;
	
	for (int i = 0; i < nc - 2; ++i){
		file >> line;
	}
	
	double **M = mat_alloc(nr,nc);
	
	for (int i = 0; i < nr; ++i){
        for (int j = 0; j < nc; ++j){
			file >> M[i][j];
		}
    }
	
	file.close();
	return M;
}
double **import_mat(std::string filename, long int &nr, int &nc){
	double line;
	std::ifstream file;
	file.open(filename);
	file >> nr >> nc;
	
	for (int i = 0; i < nc - 2; ++i){
		file >> line;
	}
	
	double **M = mat_alloc(nr,nc);
	
	for (int i = 0; i < nr; ++i){
        for (int j = 0; j < nc; ++j){
			file >> M[i][j];
		}
    }
	
	file.close();
	return M;
}
double **import_mat(std::string filename, long int &nr, long int &nc){
	double line;
	std::ifstream file;
	file.open(filename);
	file >> nr >> nc;
	
	for (int i = 0; i < nc - 2; ++i){
		file >> line;
	}
	
	double **M = mat_alloc(nr,nc);
	
	for (int i = 0; i < nr; ++i){
        for (int j = 0; j < nc; ++j){
			file >> M[i][j];
		}
    }
	
	file.close();
	return M;
}
double **import_mat(std::string filename, int &nr, long int &nc){
	double line;
	std::ifstream file;
	file.open(filename);
	file >> nr >> nc;
	
	for (int i = 0; i < nc - 2; ++i){
		file >> line;
	}
	
	double **M = mat_alloc(nr,nc);
	
	for (int i = 0; i < nr; ++i){
        for (int j = 0; j < nc; ++j){
			file >> M[i][j];
		}
    }
	
	file.close();
	return M;
}


#endif