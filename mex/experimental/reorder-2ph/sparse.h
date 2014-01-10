/*------------------------------------------------------------
 File       : sparse.h
 Created    : Thu Aug  7 16:10:06 2008
 Author     : Jostein R. Natvig <Jostein.R.Natvig@sintef.no>
 Revision   : $Id: sparse.h 136 2008-08-26 13:07:00Z jrn $


 Description


 Changes

------------------------------------------------------------*/

#ifndef _SPARSE_H_
#define _SPARSE_H_

typedef struct
{
  int     m;        /* Number of rows                                 */
  int     n;        /* Number of columns                              */
  int     nres;     /* Reserved space for non-zeros                   */
  int    *ia;       /* Position of row i                              */
  int    *ja;       /* Neighbour vertices (column number, j)          */
  double * a;       /* Flux to neighbour vertices (value in (i,j))    */
} sparse_t;

sparse_t * sparse_alloc    (int nrows, int ncols, int nnz);
void       sparse_free     (sparse_t *A);

void       sparse_permute  (sparse_t *A, int *rowP);

int        sparse_is_triangular(const sparse_t *A, const int *P);

#endif
