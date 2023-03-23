#include <stdlib.h>
#include <stdio.h>
#include <omp.h>

#include <mex.h>
#include "jacobi_basis_faster.h"
/* ---------------------------------------------------------------------- */
void
mexFunction(int nlhs,       mxArray *plhs[],
            int nrhs, const mxArray *prhs[])
/* ---------------------------------------------------------------------- */
{
    double *I;
    
    int *nI, *nC, *nF, *nO, i;
    
    double *I0, *D, *vals;
    
    int          *ii, *iiPos,
                 *nfine, *ncoarse, 
                 *bnd, *bndPos, *interact, 
                 *interactPos, *edges, *edgesPos, *overlap;
    
    I0 = mxGetPr(prhs[0]);
    nI = (int*) mxGetData(prhs[1]);
    nF = (int*) mxGetData(prhs[2]);
    nC = (int*) mxGetData(prhs[3]);
    nO = (int*) mxGetData(prhs[15]);
    
    D     = (double*) mxGetData(prhs[4]), /* Diagonal */
    ii    = (int*)    mxGetData(prhs[5]), /* ii */
    iiPos = (int*)    mxGetData(prhs[6]), /* index pos in ii */
    vals  = (double*) mxGetData(prhs[7]), /* matrix values */

    bnd    = (int*)    mxGetData(prhs[8]), /* Boundary */
    bndPos = (int*)    mxGetData(prhs[9]), /* Boundary pos */
    interact = (int*)    mxGetData(prhs[10]), /* interaction */
    interactPos = (int*)    mxGetData(prhs[11]), /* interaction pos */

    edges    = (int*)    mxGetData(prhs[12]), /* edge intersection */
    edgesPos = (int*)    mxGetData(prhs[13]), /* edge intersection pos */

    overlap = (int*)    mxGetData(prhs[14]), /* Overlap */
    
    printf("%d elements \n", nI[0]);
    /* Make output interpolator via copy (maybe excessive) */
    plhs[0] = mxCreateDoubleMatrix((mwSize)nI[0], (mwSize)1, mxREAL);
    I = mxGetPr(plhs[0]);
    
    for(i=0; i<nI[0]; i++){
        I[i] = I0[i];
    }

    mexPrintf("Hello World!\n I got %d outputs and %d inputs!\n", nlhs, nrhs);
    /*
    #pragma omp parallel
    printf("Hello from thread %d, nthreads %d\n", omp_get_thread_num(), omp_get_num_threads());
    */
    make_jacobi_basis(I, nI[0], /* Initial interpolator and size*/
                      nF[0], nC[0], /* Fine and coarse sizes    */ 
                      D, ii, iiPos, vals, bnd, bndPos, interact, interactPos,
                      edges, edgesPos, overlap,
                      nO[0] /* Overlap count */
            );
}
