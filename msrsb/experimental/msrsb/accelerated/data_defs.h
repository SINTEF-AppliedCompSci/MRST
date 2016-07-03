#pragma once
#ifndef __restrict__
#define __restrict__  
#endif
struct Grid
{
	int n_fine;
	int n_coarse;
	int n_basis;
	int * support; // n_basis elements
	int * celltypes; // n_basis elements
	int * offsets; // n_coarse + 1 elements
};

struct ConnMatrix
{
	int n_i;
	int n_conn;
	double * conn; // n_i x n_c
	int * j_index; // n_i x n_c

	int *  loc_index;  // n_basis x n_c
	double *  loc_conn; // n_basis x n_c
};
