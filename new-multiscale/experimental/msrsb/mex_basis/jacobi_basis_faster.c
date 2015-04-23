/* Written by Olav Moyner - olav.moyner@sintef.no */

#include <stdlib.h>
#include <stdio.h>

#include <math.h>
#include <mex.h>
#include <time.h>
#include <assert.h>
#include <omp.h>
#include "basis_defs.h"
#include "jacobi_basis_faster.h"

#ifndef NAN
static const unsigned long __nan[2] = {0xffffffff, 0x7fffffff};
#define NAN (*(const float *) __nan)
#endif

/* Bool indicating if selective update is enabled */
#define SELECTIVE 1
/* Default parallel chunk size */
#define LOOPCHUNK 100
/* Check convergence every CONVITS iterations */
#define CONVITS 10

/* Input arguments:
 * I            - Pointer to the memory where the prolongation operator is located.
 *
 * nf           - Number of fine cells in grid
 *
 * nc           - Number of coarse blocks in grid
 *
 * nel          - The length of the I array.
 *
 * ngb          - Total number of fine cells defined as boundary of some block
 *
 * cells        - Pointer to cell mapping for I. Same length, and cells[i]
 *               contains the fine index that will be mapped to basis value 
 *               I[i].
 *
 * bndIndicator - Array of length nel with indicator values:
 *                     0 indicates cell is not part of any boundary
 *                     1 indicates cell is locally a boundary
 *                     2 indicates cell is not locally boundary, but it is 
 *                       the same fine cell as some other coarse block's 
 *                       boundary. (Used for convergence screening)
 *
 * cellPos      - Index into cells and bndIndicator. Cells of coarse block 
 *                i can be found between from cellPos[i] to cellPos[i+1]-1.
 *
 * GB_globalCells - ngb long array indicating the global fine cell indices 
 *                  of each boundary cell.
 *
 * GBLocal        - array indicating the local indices (i.e. into 
 *                  the I matrix).
 *
 * GBLocalMap     - Map into GBLocal. For boundary element i, the local 
 *                  indices of elements in I that correspond to the global 
 *                  cell index GB_globalCells(i) is found at indices 
 *                  j between GBLocalMap(i) and GBLocalMap(i+1)-j in GBLocal.
 *
 * subsys         - Pointer array of length nc of all the local system matrices.
 *
 * tol            - Convergence tolerance (interpreted as infty norm for
 *                  all cells that are not part of any boundary.
 *
 * maxiter        - Maximum number of iterations.
 *
 * omega          - Relaxation parameter for jacobi iterations. 
 *                  Typically 1 or 2/3.
 * 
 *                  
 * 
 */
int compute_jacobi_basis(
        double *__restrict I, 
        int nf, int nc, int nel, int ngb,
        int *__restrict cells,
        int *__restrict bndIndicator, 
        int *__restrict cellPos,
        int *__restrict GB_globalCells,
        int *__restrict GBLocal, 
        int *__restrict GBLocalMap,
        int *__restrict cellIsActive,
        struct submatrix *__restrict subsys, 
        double tol, int maxiter, double omega)
{
    /* Loop indices/loop helpers */
    int i, j, k, d, pos, row, col;
    int offset, globCell, converged;

    int itcount = 0;
    
    int checkConvergence;
    
    int normalize = 1;
    
    double *blockA, *blockB, *sumval, *res, *resBlock;
    double *I_prev, *I_curr;
    
    
    double cstrength, diagstrength, newv;
    
    clock_t t0 = clock();
    clock_t elapsedtime;
    
    int msec;
    
    double rloc;
    
    mexPrintf("\nTol: %g, maxiter: %d, omega: %g\n", tol, maxiter, omega);

    /* Two blocks of memory for the nonzeros of the basis function */
    blockA  = (double*)malloc(sizeof(double)*nel);
    blockB =  (double*)malloc(sizeof(double)*nel);
    /* Residual values for each basis cell */
    resBlock = (double*)malloc(sizeof(double)*nel);
    
    /* Storage for the sum of all values used for normalization. 
     * Global size even though it is a small subset in practice 
     * to avoid index mess. */
    sumval =  (double*)malloc(sizeof(double)*nf);
    res    =  (double*)malloc(sizeof(double)*nc);

    /* First pointer set to initial basis function */
    I_prev = I;
    I_curr = blockA;
    
    /* Loop over max iterations */
    for(i = 0; i < maxiter; i++){
        /* Convergence test is somewhat expensive, only do it every now and 
         * then and only if residual increment makes any sense */
        checkConvergence = (((i+1) % CONVITS) == 0) && tol < 1;
        
        /* Re-zero memory and set current iterator to be previous values */
        #pragma omp parallel for schedule(dynamic, LOOPCHUNK)
        for(j = 0; j<nel; j++){
            I_curr[j] = I_prev[j];
            resBlock[j] = 0;
        }
        /* Do the jacobi iterations for each coarse block */
        #pragma omp parallel for private(j, k, pos, offset, d, row, col, cstrength, diagstrength, rloc)
        for(j = 0; j < nc; j++){
            offset = cellPos[j];
            d = cellPos[j+1] - cellPos[j];
            for(col = 0; col < d; col++){
                for(pos = subsys[j].rowPos[col]; pos < subsys[j].rowPos[col+1]; pos++){
                    /* Local row index */
                    row = subsys[j].row[pos];
                    
                    if(SELECTIVE){
                        if(!cellIsActive[cells[offset + row]]){
                            continue;
                        }
                    }
                    
                    if(!checkConvergence){
                        /* Compute it inline for speed if we do are not checking convergence this iteration
                         * subsys[j].values[pos] = connection strength
                         * subsys[j].diagonal[row] = diagonal entry
                         */
                        I_curr[offset + row] -= omega*subsys[j].values[pos]*I_prev[offset + col]/subsys[j].diagonal[row];
                    }else{
                        /* Compute residual in resBlock, which is reduced 
                         * to each coarse block later on */
                        rloc = subsys[j].values[pos]*I_prev[offset + col]/subsys[j].diagonal[row];
                        resBlock[offset + row] += rloc;
                        /* Set actual basis value by relaxing with omega factor */
                        I_curr[offset + row] -= omega*rloc;
                    }
                }
            }
            /* Inf norm check over all cells not part of some boundary */
            if(checkConvergence){
                res[j] = 0;
                for(k = 0; k < d; k++){
                    if(bndIndicator[k + offset] < 1){
                        res[j] = max(res[j], fabs(resBlock[offset + k]));
                    }
                }
            }
        }
        
        /* Normalize */
        if(normalize){
            #pragma omp parallel for schedule(dynamic, LOOPCHUNK)
            for(j = 0; j < nf; j++){
                /* Set memory to zero */
                sumval[j] = 0;
            }
            #pragma omp parallel for schedule(dynamic, LOOPCHUNK)
            for(j = 0; j < nel; j++){
                /* Flag 1 indicates boundary of local cell, normalize 
                 * and add to removed values */
                
                if(SELECTIVE){
                    if(!cellIsActive[cells[j]]){
                        continue;
                    }
                }
                if(bndIndicator[j] == 1){
                    #pragma omp atomic
                    sumval[cells[j]] += I_curr[j];
                    I_curr[j] = 0;
                }
            }
            #pragma omp parallel for private(globCell, k)
            for(j = 0; j < ngb; j++){
                globCell = GB_globalCells[j];
                
                if(SELECTIVE){
                    if(!cellIsActive[globCell]){
                        continue;
                    }
                }
                for(k = GBLocalMap[j]; k < GBLocalMap[j+1]; k++){
                    /* Re-normalize all the values based on the removed
                     * boundary values (sum is 1 previous to removal of 
                     * values, so divide by 1- (removed values)) */
                    newv = I_curr[GBLocal[k]]/(1 - sumval[globCell]);
                    I_curr[GBLocal[k]] = min(max(newv, 0), 1);
                }
            } 
        } 
        
        /* Memory management - switch previous/current pointers around */
        if(i % 2 == 0){
            /* Write to blockA on odd passes */
            I_prev = blockA;
            I_curr = blockB;
        }else{
            I_prev = blockB;
            I_curr = blockA;
        }
        itcount = i + 1;
        
        if(checkConvergence){
            /* We converge globally if all coarse blocks have 
             * converged locally*/
            converged = 1;
            for(j = 0; j < nc; j++){
                converged = converged && (res[j] <= tol);
            }
            if(converged){
                break;
            }
        }
        
        /* Renormalize completely every 100 iterations or so just in case 
         * of drifting in the range of values. Typically this happens if we
         * iterate more than 1000 iterations (but safeguard against it anyway).
         */
        if(i % 100 == 0){
            renormalize(sumval, I_prev, cells, nel, nf);
        }
        
    }
    /* A GLOBAL re-normalization to take care of any rounding errors */
    renormalize(sumval, I_prev, cells, nel, nf);
    
    /* Write to the external pointer and return the result */
    #pragma omp parallel for
    for(i = 0; i< nel; i++){
        I[i] = I_prev[i];
    }

    
    /* Free allocated memory */
    free(blockA);
    free(blockB);
    free(sumval);
    free(res);
    free(resBlock);
    
    elapsedtime = clock() - t0;
    msec = elapsedtime * 1000 / CLOCKS_PER_SEC;
    
    printf("Computed result with %d iterations in CPU time %d seconds %d milliseconds (%d ms per iteration)\n",
            itcount,  msec/1000, msec%1000, msec/max(itcount, 1));
    return 1;
}

int renormalize(double *sumval, double *I_curr, int *cells, int nel, int nf){
    int j;
    #pragma omp parallel for schedule(dynamic, LOOPCHUNK)
    for(j = 0; j < nf; j++){
        /* Set memory to zero */
        sumval[j] = 0;
    }
    #pragma omp parallel for schedule(dynamic, LOOPCHUNK)
    for(j = 0; j < nel; j++){
        #pragma omp atomic
        sumval[cells[j]] += I_curr[j];
    }
    #pragma omp parallel for schedule(dynamic, LOOPCHUNK)
    for(j = 0; j < nel; j++){
        I_curr[j] = I_curr[j]/sumval[cells[j]];
    }
    return 1;
}

