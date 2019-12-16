/*
 * Copyright 2016 (c) SINTEF ICT, Applied Mathematics.
 * Maintained by Olav Moyner <olav.moyner@sintef.no>
 */
#include <iostream>
#include "basis_solver.h"
#include <string>
#include <chrono>
#include <math.h>
#include <algorithm>

#include <omp.h>
#include "mex_basis_solver.h"

#ifdef USEMEX
#include <mex.h>
#include <matrix.h>
#endif

#include <cstring>
#define TIME_NOW std::chrono::high_resolution_clock::now();

int main(int argc, char* argv[])
{
	if(argc != 2){
		std::cout << "Please provide a single directory in standalone mode. You provided "<< argc << "\n";
		return 1;
	}
	int retval = flatFileBasis(argv[1]);
	if(retval != 0){
		std::cout << "Unable to read necessary files. Aborting.\n";
	}
};


int flatFileBasis(std::string fn) {
	Grid grid;
	ConnMatrix mat;

	std::string infile, outfile;
	#ifdef _WIN32
		infile = fn + "\\input\\";
		outfile = fn + "\\output\\";
	#else
		infile = fn + "/input/";
		outfile = fn + "/output/";
	#endif
	if(readInfo(&grid, infile)){return 1;};
	if(readConnMatrix(&grid, &mat, infile)){return 1;};

	// Read in basis operator
	double * basis;
	basis = new double[grid.n_basis]();
	if(readBasisOperator(&grid, basis, infile)){return 1;};

	mat.loc_index = new int[grid.n_basis*mat.n_conn]();
	mat.loc_conn = new double[grid.n_basis*mat.n_conn]();

	getBasis(&grid, &mat, basis, -1, 100, 0.66);
	// Write output
	writeBasisOperator(&grid, outfile, basis);

	delete[] grid.support;
	delete[] grid.celltypes;
	delete[] grid.offsets;

	delete[] basis;

	delete[] mat.conn;
	delete[] mat.j_index;

	delete[] mat.loc_index;
	delete[] mat.loc_conn;
	return 0;
};
/* ---------------------------------------------------------------------- */

#ifdef USEMEX
void mexFunction(
	int          nlhs,
	mxArray      *plhs[],
	int          nrhs,
	const mxArray *prhs[]
	)
{
	/* If we have recieved an empty array, fall back to flatfiles */
	printf("Recieved %d input arguments and requested %d outputs\n", nrhs, nlhs);

	if (nrhs == 1) {
		int retval = flatFileBasis(mxArrayToString(prhs[0]));
		if(retval){
			printf("Unable to read necessary files. Aborting.\n");
			return;
		}
		printf("Flatfile processing complete. Basis written to file...\n");
		return;
	}
	if(nrhs != 9){
		return;
	}
	Grid grid;
	ConnMatrix mat;
	double * basis;
	std::string fn;


	// offsets, support, types, mat, j_index, I
	/* Offsets for each coarse block */
	grid.offsets = (int*)mxGetData(prhs[0]);

	/* Number of coarse blocks */
	const size_t * offset_len = mxGetDimensions(prhs[0]);
	grid.n_coarse = offset_len[0] - 1;

	/* Number of values in basis function */
	grid.n_basis = grid.offsets[grid.n_coarse];

	/* n_el long array, mapping each basis cell into fine cells */
	grid.support = (int*)mxGetData(prhs[1]);

	/* Types of each element in the support */
	grid.celltypes = (int*)mxGetData(prhs[2]);


	/* Total number of fine cells */
	const size_t * matrix_dim = mxGetDimensions(prhs[3]);
	grid.n_fine = matrix_dim[1];

	/* Building matrix */
	auto t1 = TIME_NOW;
	buildMatrixFromMxSparse(&mat, prhs[3]);
	auto t2 = TIME_NOW;
	std::cout << "Matlab matrix processed in "
		<< std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count()
		<< " ms"  << std::endl;

	/* Initial basis guess */
	double * initial_basis = (double*)mxGetData(prhs[4]);

	double tol = mxGetScalar(prhs[5]);
	int maxiter = mxGetScalar(prhs[6]);
	double omega = mxGetScalar(prhs[7]);
    int given_threads = mxGetScalar(prhs[8]);
    int max_threads = std::min(given_threads, omp_get_num_procs());

    
	plhs[0] = mxCreateDoubleMatrix((mwSize)grid.n_basis, (mwSize)1, mxREAL);
	basis = (double *)mxGetData(plhs[0]);

	for (int i = 0; i < grid.n_basis; i++) {
		basis[i] = initial_basis[i];
	}
	printf("Discovered problem with %d cells and %d blocks\n", grid.n_fine, grid.n_coarse);
	printf("Problem has %d equations with up to %d connections and support in %d cells total.\n", mat.n_i, mat.n_conn, grid.n_basis);
    printf("Max number of threads: %d\n", max_threads);
    omp_set_num_threads(max_threads);
	
	mat.loc_index = new int[grid.n_basis*mat.n_conn];
	mat.loc_conn = new double[grid.n_basis*mat.n_conn];

	getBasis(&grid, &mat, basis, tol, maxiter, omega);

	delete[] mat.conn;
	delete[] mat.j_index;
	delete[] mat.loc_index;
	delete[] mat.loc_conn;
}

void buildMatrixFromMxSparse(ConnMatrix * mat, const mxArray * A) {
	//int n_el = mxGetNzmax(A);
	int n, m;
	n = mxGetN(A);
	m = mxGetM(A);
	double * data = (double *)mxGetPr(A);
	mwIndex * irs = (mwIndex *) mxGetIr(A);
	mwIndex * jcs = (mwIndex *) mxGetJc(A);

	int n_el = 0;
	for (int i = 0; i < m; i++) {
		n_el = std::max(n_el, (int)(jcs[i+1] - jcs[i]) - 1);
	}
	// printf("Matrix is %d by %d with maximum elements %d and %d nz\n", n, m, n_el, jcs[m]);

	int counts = n_el*n;
	mat->conn = new double[counts];
	mat->j_index = new int[counts];
	for (int i = 0; i < counts; i++) {
		mat->conn[i] = 0;
		mat->j_index[i] = -1;
	}

	double * diag = new double[m];
	// Massage matrix
	for (int i = 0; i < m; i++) {
		int ctr = 0;
		for (int j = jcs[i]; j < jcs[i + 1]; j++) {
			if (i == irs[j]) {
				diag[i] = data[j];
				continue;
			}
			mat->conn[i*n_el + ctr]    = data[j];
			mat->j_index[i*n_el + ctr] = irs[j];
		    ctr++;
		}
		
	}
	// Scale by missing diagonal
    #pragma omp parallel for
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n_el; j++) {
			mat->conn[i*n_el + j] = mat->conn[i*n_el + j] / diag[i];
			//printf("%f\t", mat->conn[i*n_el + j]);
		}
		//printf("\n");
	}

	mat->n_conn = n_el;
	mat->n_i = n;

	delete[] diag;
	return;
}
#endif
void printMatrix(ConnMatrix * mat){
	printf("J_indices:\n");
	for (int i = 0; i < mat->n_i; i++) {
		for(int j = 0; j < mat->n_conn; j++){
			printf("%d ", mat->j_index[(mat->n_conn)*i + j]);
		}
		printf("\n");
	}
	printf("Connections:\n");
	for (int i = 0; i < mat->n_i; i++) {
		for(int j = 0; j < mat->n_conn; j++){
			printf("%1.2g ", mat->conn[(mat->n_conn)*i + j]);
		}
		printf("\n");
	}
}

void printGrid(Grid * grid){
	int nc = grid->n_coarse;
	int nb = grid->n_basis;

	printf("Support:\n");
	for (int i = 0; i < nb; i++) {
		printf("%d\n", grid->support[i]);
	}

	printf("Types:\n");
	for (int i = 0; i < nb; i++) {
		printf("%d\n", grid->celltypes[i]);
	}

	printf("Offsets:\n");
	for (int i = 0; i < nc + 1; i++) {
		printf("%d\n", grid->offsets[i]);
	}
}

void printBasis(Grid * grid, double * basis){
	printf("Basis:\n");
	for (int i = 0; i < grid->n_basis; i++) {
		printf("%1.2g\n", basis[i]);
	}
}
