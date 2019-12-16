/*
 * Copyright 2016 (c) SINTEF ICT, Applied Mathematics.
 * Maintained by Olav Moyner <olav.moyner@sintef.no>
 */
#include <iostream>
#include "basis_solver.h"
#include <fstream>
#include <iomanip>
#include <omp.h>
#include <chrono>
#include <math.h>
#include <algorithm>

#include <vector>
// Chunk size for paralleization
#define LOOPCHUNK 100
// How often we check convergence
#define CHECKITS 25
// Renormalize globally every now and then due to rounding errors (only 
// applicable for large number of iterations)
#define NORM_ITER 100
// Macro for time taking
#define TIME_NOW std::chrono::high_resolution_clock::now();
#define _OPENMP_NOFORCE_MANIFEST 
#define INTERACTION_SORTED 1
int i, j, k, l;

void getBasis(Grid * grid, ConnMatrix * mat, double * basis, double tol, int N, double relax) {
	/* Set up matrix operators */
	std::cout << "Creating mappings..." << std::endl;
	auto t1 = TIME_NOW;

	setupCoarseMapping(grid, mat);

	std::cout << "Starting basis generation..." << std::endl;
	auto t2 = TIME_NOW;

	/* Generate basis functions */
	computeBasis(grid, mat, basis, tol, N, relax);


	auto t3 = TIME_NOW;
	std::cout << "Basis generation done in "
		<< std::chrono::duration_cast<std::chrono::milliseconds>(t3 - t1).count()
		<< " ms"  << std::endl;
	std::cout << "---- Mapping generation: " << std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count() << "ms" << std::endl;;
	std::cout << "---- Basis functions   : " << std::chrono::duration_cast<std::chrono::milliseconds>(t3 - t2).count() << "ms" << std::endl;;
	#pragma omp barrier
}

int computeBasis(Grid * grid, ConnMatrix * mat, double * __restrict__ basis, double tol, int maxiter, double omega) {
	// Number of non-zero elements in the basis (the values explicitly assigned in sparse matrix
    int n_el = grid->n_basis;
    // Number of coarse blocks
	int n_c = grid->n_coarse;
    // Number of fine cells in total
	int n_f = grid->n_fine;
    // Upper bound on the number of off diagonal entries in a matrix row
	int n_conn = mat->n_conn;

    // Used for sum of basis functions for renormalization
	double * update_sum = new double[n_f]();
    // Storage for the update values
	double * basis_update = new double[n_el]();

	std::cout << "Starting smoothing iterations... ";
	// Iterate up to max iterations
	int iter = 0;
	
	while(iter < maxiter){
        // We check convergence every CHECKITS and as long as we have not converged.
		bool check_convergence = (iter % CHECKITS) == 0 && iter > 0 && tol > 0;
		bool done = check_convergence;
		// Iterate over each coarse block in parallel
		#pragma omp barrier
		#pragma omp parallel for
		for (int coarseIx = 0; coarseIx < n_c; coarseIx++) {
            // Compute the position in the sparse matrix where we have cells
			int start = grid->offsets[coarseIx];
			int stop  = grid->offsets[coarseIx+1];
			// Iterate over local support
			for (int j = start; j < stop; j++) {
				basis_update[j] = basis[j];
				// Iterate over local values
				for (int k = 0; k < n_conn; k++) {
					int jx = mat->loc_index[j*n_conn + k];
					if (jx == -1) {
                        // We have reached the end of the active connections
						break;
					}
                    // This is the tentative update for each basis element
					basis_update[j] +=  mat->loc_conn[j*n_conn + k]* basis[jx];
				}
			}
		}
        // Modify the update based on type of cell (inner, boundary or global boundary)
		#pragma omp barrier
		#pragma omp parallel for
		for (int i = 0; i < n_el; i++) {
			int type = grid->celltypes[i];
			if(type == 0){
				// The cell is category 1: We simply update the value.
				basis[i] -= omega*basis_update[i];
                if (check_convergence && fabs(basis_update[i]) > tol){
                    done = false;
                }
            }
			else if (type == 2) {
                // The cell is part of the boundary for some other cell. We
                // need to keep track of what the nonzero values are.
				basis[i] -= omega*basis_update[i];
                #pragma omp atomic
                update_sum[grid->support[i]] += basis[i];
			}
		}
		#pragma omp barrier
		#pragma omp parallel for
		for (int i = 0; i < n_el; i++) {
            // Fix partition of unity for cells that are part of some boundary.
			if (grid->celltypes[i] == 2) {
				double us = update_sum[grid->support[i]];
				basis[i] = basis[i]/us;
			}
		}
		#pragma omp barrier
		#pragma omp parallel for
		for (int i = 0; i < n_f; i++) {
			update_sum[i] = 0;
		}
		#pragma omp barrier
		if (done) {
			break;
		}
		iter++;
	}

	delete[] update_sum;
	delete[] basis_update;
	std::cout << "done after " << iter << " iterations." << std::endl;
	return 0;
}

void renormalize(Grid * grid, ConnMatrix * mat, double * basis){
    // Renormalize so that the basis functions have sum one globally
	int n_el = grid->n_basis;
	int n_f = grid->n_fine;
	double * sumval = new double[n_f]();
	#pragma omp parallel for
	for(int i = 0; i < n_el; i++){
		#pragma omp atomic
		sumval[grid->support[i]] += basis[i];
	}
	#pragma omp barrier
	#pragma omp parallel for
	for(int i = 0; i < n_el; i++){
			basis[i] = basis[i]/sumval[grid->support[i]];
	}
	#pragma omp barrier
	delete[] sumval;
}

int readInfo(Grid * grid, std::string file_path) {
	std::ifstream infofile, supportfile, typefile;
	if(openFlatFile(file_path + "info.txt", &infofile)){return 1;};
	if(openFlatFile(file_path + "support.txt", &supportfile)){return 1;};
	if(openFlatFile(file_path + "types.txt", &typefile)){return 1;};

	int nf;
	int nc;
	int n_basis;

	// Read sizes
	infofile >> nf;
	infofile >> nc;
	// Read offsets
	grid->offsets = new int[nc + 1];
	int i = 0;
	while (!infofile.eof() && i < nc + 1) {
		infofile >> grid->offsets[i];
		i++;
	}
	
	n_basis = grid->offsets[nc];
	// Read support
	grid->support = new int[n_basis];
	grid->celltypes = new int[n_basis];
	i = 0;
	while (!supportfile.eof() && !typefile.eof() && i < n_basis) {
		supportfile >> grid->support[i];
		//std::cout << i << " -> " << grid->support[i] << "\r\n";
		typefile >> grid->celltypes[i];
		i++;
	}
	if (i != n_basis) {
		std::cout << "Something went wrong\n";
		std::cout << "Got " << i << " should be " << n_basis;
		return 1;
	}

	// store misc constants
	grid->n_basis = n_basis;
	grid->n_fine = nf;
	grid->n_coarse = nc;

	return 0;
};

int readConnMatrix(Grid * grid, ConnMatrix * mat, std::string file_path) {
	std::ifstream f, fs;
	if(openFlatFile(file_path + "matrix.txt", &f)){return 1;};
	if(openFlatFile(file_path + "sparsity.txt", &fs)){return 1;};
	int n_conn, n_el;
	n_el = grid->n_fine;
	f >> n_conn;



	mat->conn    = new double[n_el*n_conn];
	mat->j_index = new int[n_el*n_conn];
	mat->n_i = n_el;
	mat->n_conn = n_conn;
	std::cout << "Matrix is " << mat->n_i << " by " << n_conn << "\r\n";

	int i = 0;
	int ix;
	while (!f.eof() && i < n_el) {
		for (int j = 0; j < n_conn; j++) {
			ix = i*n_conn + j;
			f >> (mat->conn)[ix];
			fs >> (mat->j_index)[ix];
		}
		i++;
	}

	if (0) {
		for (int i = 0; i < n_el; i++) {
			for (int j = 0; j < n_conn; j++) {
				ix = i*n_conn + j;
				std::cout << mat->conn[ix] << " ";
			}
			std::cout << "\r\n";
		}
	}
	return 0;
};

int readBasisOperator(Grid * grid, double * basis, std::string file_path) {
	int n_el = grid->n_basis;

	std::ifstream f;
	openFlatFile(file_path + "operator.txt", &f);

	int i = 0;
	while (!f.eof() && i < n_el) {
		f >> basis[i];
		i++;
	}
	return 0;
};

void setupCoarseMapping(Grid * grid, ConnMatrix * mat){
	int n_el = grid->n_basis;
	int n_c =  grid->n_coarse;
	int n_conn = mat->n_conn;

	for (int i = 0; i < n_el*n_conn; i++) {
		mat->loc_index[i] = -1;
		mat->loc_conn[i] = 0;
	}
	#pragma omp parallel for
	for (int coarseIx = 0; coarseIx < n_c; coarseIx++) {
		int start = grid->offsets[coarseIx];
		int stop  = grid->offsets[coarseIx + 1];
		for (int j = start; j < stop; j++) {
			int c = grid->support[j];
			for (int k = 0; k < n_conn; k++) {
				int c_n = mat->j_index[c*n_conn + k];
				int l;
				bool done = false;
				if(c_n == -1){
					break;
				}
				if(INTERACTION_SORTED){
					// This branch assumes support is sorted.
					if (c_n > c) {
						for (l = j; l < stop; l++) {
							if (grid->support[l] == c_n) {
								done = true;
								break;
							}
						}
					}
					else {
						for (l = j - 1; l >= start; l--) {
							if (grid->support[l] == c_n) {
								done = true;
								break;
							}
						}
					}
				}else{
					for (l = start; l < stop; l++) {
						if (grid->support[l] == c_n) {
							done = true;
							break;
						}
					}
				}
                /*
				if(!done && grid->celltypes[j] != 1){
					printf("Error! %d: %d (mapping cell %d to %d) with category %d...\n", coarseIx, j, c, c_n, grid->celltypes[j]);
				}
                */
				if(done){
					mat->loc_index[j*n_conn + k] = l;
					mat->loc_conn[j*n_conn + k] = mat->conn[c*n_conn + k];
				}
			}
		}
	}
#pragma omp barrier
}

int writeBasisOperator(Grid * grid, std::string file_path, double * basis) {
	int n_el = grid->n_basis;

	std::ofstream myfile(file_path + "operator.txt");

	if (myfile.is_open())
	{
		for (int i = 0; i < n_el; i++) {
			myfile << std::fixed << std::setprecision(8) << basis[i];
			myfile << " ";
		}
		return 0;
	}else{
		return -1;
	}
};

int openFlatFile(std::string file_path, std::ifstream * file) {
	file->open(file_path.c_str());
	if (!(file->is_open())) {
		std::cout << "Unable to open file: " << file_path << "\r\n";
		return 1;
	}
	return 0;
}

