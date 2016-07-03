#include <stdlib.h>
#include <stdio.h>

#include <math.h>
#include <mex.h> 
#include <time.h>

#include <omp.h>

#include "jacobi_basis_faster.h"

#ifndef NAN
    static const unsigned long __nan[2] = {0xffffffff, 0x7fffffff};
    #define NAN (*(const float *) __nan)
#endif


 int make_jacobi_basis(double* I, int nI, int nF, int nC, 
                       double* diag, int* ii, int* ipos, double* vals, 
                       int* bnd, int* bndpos, int* interact, int* interactpos,
                       int* edge, int* edgepos,
                       int* overlap, int nOverlap)
{
    int i, j, k, itcount, nel, offset, w, ci, cn, locind, k1, k2, ip, b, ok;
    int maxiter = 100;
    
    int normalize = 1;
    
    double val, d;
    
    double *blockA, *blockB, *sumval;
    double *I_prev, *I_curr;
    
    double *res, tmp;
    
    double tol  = 1e-2;
    double omega = 0.666;
    
    clock_t t0 = clock();
    clock_t time;
    
    int msec;
    
    blockA  = (double*)malloc(sizeof(double)*nI);
    blockB =  (double*)malloc(sizeof(double)*nI);
    
    sumval =  (double*)malloc(sizeof(double)*nF);
    res    = (double*)malloc(sizeof(double)*nC);
    
    I_prev = I;
    
    I_curr = blockA;
    
    for(i = 0; i<nI; i++){
        I_curr[i] = 0;
    }
    itcount = maxiter;
    /* Loop over max iterations */
    for(i = 0; i<maxiter; i++){
        /* #pragma omp parallel for schedule(dynamic, 10) */ 
        for(j = 0; j < nC; j++)
        {
            if(i == 0 && j == 0){
                printf("Solving iterations using %d threads.\n", omp_get_num_threads());
            }
            res[j] = 0;
            /* Loop over all cells */
            for(k1 = interactpos[j]; k1 < interactpos[j+1]; k1++){
                ci = interact[k1];
                val = diag[ci]*I_prev[k1];
                for(ip = ipos[ci]; ip < ipos[ci+1]; ip++){
                    cn = ii[ip];
                    k2 = find_cell_neighbor(interactpos, interact, j, cn, ci, k1);
                    if(k2 < 0){
                        /* Out of bounds */
                        continue;
                    }
                    val = val + I_prev[k2]*vals[ip];
                }
                I_curr[k1] = I_prev[k1] - omega*val/diag[ci];
                res[j] = res[j] + fabs(val/diag[ci]);
            }
        }
        /* #pragma omp barrier */
        if(normalize){
            /* Store the removed values somewhere */
            for(j = 0; j < nF; j++){
                sumval[j] = 0;
            }
            /* Set zero values at each local boundary, using bookkeeping to store the removed values */
            for(j = 0; j < nC; j++){
                /* Set local boundary to zero */
                for(b = bndpos[j]; b < bndpos[j+1]; b++){
                    k2 = find_cell_neighbor_simple(interactpos, interact, j, bnd[b]);
                    sumval[bnd[b]] += I_curr[k2];
                    I_curr[k2] = 0;
                }
            }
            /* Renormalize the edge values */
            /*
             #pragma omp barrier
             #pragma omp parallel for schedule(dynamic, 10)
             */
            for(j = 0; j < nC; j++){
                for(k1 = edgepos[j]; k1 < edgepos[j+1]; k1++){
                    ci = edge[k1];
                    if(sumval[ci] > 0){
                        k2 = find_cell_neighbor_simple(interactpos, interact, j, ci);
                        if(sumval[ci] < 1){
                           /* We know that sum of all cells is 1, so subtract any removed values */
                           I_curr[k2] = I_curr[k2]/(1-sumval[edge[k1]]);
                        }
                    }
                }
            }
        }
        /* Memory management - put things in place */
        if(i % 2 == 0){
            /* Write to blockA on odd passes */
            I_prev = blockA;
            I_curr = blockB;
        }else{
            I_prev = blockB;
            I_curr = blockA;
        }
        
        /* Commented out - need to handle residual near edges in a nice manner */
        /*
        ok = 1;
        tmp = 0;
        for(j = 0; j < nC; j++)
        {
            nel = (interactpos[j+1] - interactpos[j]);
            ok = ok && (res[j]/(double)nel < tol);
            tmp = max(tmp, res[j]/((double)nel));
        }
        mexPrintf("Biggest value: %g\n", tmp);
        if(ok){
            itcount = i + 1;
            break;
        } */
    }
    for(i = 0; i<nI; i++){
        I[i] = I_prev[i];
    }
    free(blockA);
    free(blockB);
    free(sumval);
    free(res);
    
    time = clock() - t0;
    msec = time * 1000 / CLOCKS_PER_SEC;
    printf("Computed result with %d iterations in %d seconds %d milliseconds\n", itcount,  msec/1000, msec%1000);
    return 1;
}
 
int find_cell_neighbor(int *interactpos, int *interact, int coarse, int target, int current, int currentk){
    int k;
    /* Seek upwards in interaction */
    if(target > current){
        for(k = currentk+1; k < interactpos[coarse+1]; k++){
            if(interact[k] == target){
                return k;
            }
        }
    }
    else /*Seek downwards in interaction */
    {
        for(k = currentk-1; k > interactpos[coarse]-1; k--){
            if(interact[k] == target){
                return k;
            }
        }
    }
    /* We did not find the connection, return error */
    return -1;
}

int find_cell_neighbor_simple(int *interactpos, int *interact, int coarse, int target){
    int k;
    for(k = interactpos[coarse]; k < interactpos[coarse+1]; k++){
        if(interact[k] == target){
            return k;
        }
    }
}
