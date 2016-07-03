#include <stdlib.h>
#include <stdio.h>
#include <omp.h>

#include <mex.h>

#include "basis_defs.h"
#include "jacobi_basis_faster.h"

/* ---------------------------------------------------------------------- */
void
mexFunction(int nlhs,       mxArray *plhs[],
            int nrhs, const mxArray *prhs[])
/* ---------------------------------------------------------------------- */
{
    double *I, *I0;
    
    int  nf,              nc,       nel,     ngb,
        *cells,          *isBnd,   *cellPos, *cellIsActive,
        *GB_globalCells, *GBLocal, *GBLocalMap;
    
    int i;
    
    double tol, omega;
    int maxiter;

    struct submatrix *subsys; 
    /* Sizes */
    nf = (int) mxGetScalar(prhs[0]);
    nc = (int) mxGetScalar(prhs[1]);
    nel = (int) mxGetScalar(prhs[2]);
    
    cells = (int*) mxGetData(prhs[4]);
    isBnd = (int*) mxGetData(prhs[5]);
    cellPos = (int*) mxGetData(prhs[6]);
    
    GB_globalCells = (int*) mxGetData(prhs[7]);
    GBLocal        = (int*) mxGetData(prhs[8]);
    GBLocalMap     = (int*) mxGetData(prhs[9]);
    
    cellIsActive     = (int*) mxGetData(prhs[14]);
    
    tol = (double) mxGetScalar(prhs[15]);
    maxiter = (int) mxGetScalar(prhs[16]);
    omega = (double) mxGetScalar(prhs[17]);
    
    ngb = (int) mxGetM(prhs[7]);
    /* 
     * mexPrintf("Global boundary has %d elements \n", ngb);
     * mexPrintf("I got %d outputs and %d inputs!\n", nlhs, nrhs);
     * mexPrintf("Nf: %d, Nc: %d, Nel: %d\n", nf, nc, nel);
    */
    
    /* Deal with submatrices by putting them into structs */ 
    /* struct submatrix subsys[nc]; */
    subsys = (submatrix*) malloc(sizeof(submatrix)*nc);
    for(i=0; i<nc; i++){
        subsys[i].row    = (int*)    mxGetPr(mxGetCell(prhs[10], i));
        subsys[i].rowPos = (int*)    mxGetPr(mxGetCell(prhs[11], i));
        subsys[i].values = (double*) mxGetPr(mxGetCell(prhs[12], i));
        subsys[i].diagonal = (double*) mxGetPr(mxGetCell(prhs[13], i));
        subsys[i].n = (int) mxGetM(mxGetCell(prhs[10], i));
    }
    
    I0 = mxGetPr(prhs[3]);
    
    plhs[0] = mxCreateDoubleMatrix((mwSize)nel, (mwSize)1, mxREAL);
    I = mxGetPr(plhs[0]);

    for(i=0; i<nel; i++){
        I[i] = I0[i];
    }
    compute_jacobi_basis(I, nf, nc, nel, ngb, 
         cells, isBnd, cellPos, 
         GB_globalCells, GBLocal, GBLocalMap, cellIsActive, 
         subsys, tol, maxiter, omega);
}
