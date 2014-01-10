/*------------------------------------------------------------
 File       : sparse.c
 Created    : Thu Aug  7 16:10:33 2008
 Author     : Jostein R. Natvig <Jostein.R.Natvig@sintef.no>
            : Knut-Andreas Lie <Knut-Andreas.Lie@sintef.no>
 Revision   : $Id: sparse.c 186 2008-10-29 16:52:37Z jrn $

 Description: module for handling sparse matrices
------------------------------------------------------------*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include "sparse.h"

#if MATLAB_MEX_FILE
#include <mex.h>
#define print   mexPrintf
#define malloc  mxMalloc
#define calloc  mxCalloc
#define realloc mxRealloc
#define free    mxFree
#else

#endif


sparse_t *sparse_alloc(int nrows, int ncols, int nnz)
/*---------------------------------------------------------------------------
 * Input arguments:
 *
 *     nrows  - number or rows
 *     ncols  - number of columns
 *     nnz    - number of nonzero elements
 *
 * Output:
 *
 *     return - Sparse nrows x ncols matrix witgh room for nnz nonzeros.
 *---------------------------------------------------------------------------*/
{
  assert(nrows >= 0);
  assert(ncols >= 0);
  assert(nnz >= 0);
  sparse_t *A;
  if(!(A = malloc(sizeof(sparse_t))))
  {
    printf("Allocation failure in allocate_sparse - aborting.");
    exit(1);
  }

  A->m     = nrows;
  A->n     = ncols;
  A->nres  = nnz;
  if(!(A->ia = malloc((nrows+1)*sizeof(int))))
  {
    printf("Allocation failure in allocate_sparse - aborting.");
    exit(1);
  }
  /* Ensure that the size of a zero-dimensional matrix can be looked up. */
  A->ia[0] = 0;

  if(!(A->ja = malloc(nnz*sizeof(int))))
  {
    printf("Allocation failure in allocate_sparse - aborting.");
    exit(1);
  }

  if(!(A->a   = malloc(nnz*sizeof(double))))
  {
    printf("Allocation failure in allocate_sparse - aborting.");
    exit(1);
  }

  return A;
}


void sparse_free (sparse_t *A)
/*---------------------------------------------------------------------------
 * Input arguments:
 *     A  - sparse matrix allocated with alloc_sparse
 *---------------------------------------------------------------------------*/
{
  if (!A) return;

  if (A->ia) {
    free (A->ia);
    A->ia = NULL;
  }
  if (A->ja) {
    free (A->ja);
    A->ja = NULL;
  }
  if (A->a)  {
    free (A->a);
    A->a = NULL;
  }
  free(A);
  return;
}


void sparse_permute(sparse_t *A, int *P)
/*---------------------------------------------------------------------------
 * Permute rows and columsn of sparse matrix such that A_{i,j} = B_{P[i],P[j]}
 *
 * Input arguments:
 *     A  - pointer to quadratic sparse matrix
 *     P  - permutation
 *
 * Output:
 *     A  - row- and column permuted sparse matrix
 *---------------------------------------------------------------------------*/
{
  int    i,q;
  int k;
  int    m = A->m;
  if (m != A->n) {
      print("Error, cannot permute non-square matrix\n");
      exit(1);
    }

  sparse_t *B   = sparse_alloc(m, m, A->ia[m]);
  int      *invP = calloc(m, sizeof(int));

  /* Invert permutation */
  for (i=0; i<m; ++i)
    invP[P[i]] = i;

  /* Permuted row positions */
  B->ia[0]=0;
  for (i=0; i<m; ++i)
    B->ia[i+1] = B->ia[i] + A->ia[P[i]+1]-A->ia[P[i]];

  /* Permute columns */
  for (i=0; i<m; ++i) {
    const int pi = P[i];
    for (k=A->ia[pi], q=B->ia[i]; k<A->ia[pi+1]; ++k, ++q) {
      int j    = A->ja[k];
      B->ja[q] = invP[j];
      B->a [q] = A->a[k];
    }
  }

  /* Let A take over data structures of B */
  free (A->ia); A->ia = B->ia;
  free (A->ja); A->ja = B->ja;
  free (A->a ); A->a  = B->a;
  free (B);
  free (invP);
}

int sparse_is_triangular(const sparse_t *A, const int *P)
/*---------------------------------------------------------------------------
 * Check if the permuted matrix A is lower triangular
 *
 * Input arguments:
 *     A        - sparse matrix allocated with alloc_sparse
 *     P        - permutation vector
 *
 * Output:
 *     return   - 1 if A is strictly lower triangular
 *---------------------------------------------------------------------------*/
{
  int n = A->m;
  int i,k;

  if (P==NULL) {
    for (i=0; i<n; ++i)
      for (k=A->ia[i]; k<A->ia[i+1]; ++k)
	if (A->ja[k]>i)
	  return 0;
  }
  else {
    for (i=0; i<n; ++i)
      for (k=A->ia[P[i]]; k<A->ia[P[i]+1]; ++k)
	if (A->ja[k]>P[i])
	  return 0;
  }
  return 1;
}
