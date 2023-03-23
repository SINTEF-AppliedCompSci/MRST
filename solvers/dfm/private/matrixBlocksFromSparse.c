#include <stdio.h>
#include <mex.h>
#include <matrix.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <lapack.h>

/* Extract blocks on the diagonal of sparse matrix.  This is a workaround 
 * since indexing of the form A(k) does not work in Matlab for very large
 * matrices.
 */
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
  if ((nrhs != 2) || (nlhs != 1))
    {
      fprintf(stderr, "return early\n");
      return;
    }
  int i,j,k,p;

  /* Get input */
  if (!mxIsSparse(prhs[0]))
  {
     mexErrMsgTxt("First argument must be a sparse matrix.");  
  }
  
  int n = mxGetM(prhs[0]);
  if(mxGetN(prhs[0]) != n)
  {
    mexErrMsgTxt("Matrix must be square");     
  }
  
  double  *sa  = mxGetPr(prhs[0]);
  mwIndex *ia  = mxGetJc(prhs[0]);
  mwIndex *ja  = mxGetIr(prhs[0]);
  
  if(!mxIsInt32(prhs[1]))
    {
      mexErrMsgTxt("Second argument (vector of sizes) must be int32.");
    }
  int     *sz  = (int*)mxGetData(prhs[1]);
  int      nsz = mxGetNumberOfElements(prhs[1]);
   

  /* Create output vector */
  int sumsz=0;
  for (i=0; i<nsz; ++i) sumsz += sz[i];

  if (sumsz != mxGetM(prhs[0]))
  {
     mexErrMsgTxt("Block sizes does not match dimension of sparse matrix.");
  }
  
  sumsz=0;
  for (i=0; i<nsz; ++i) sumsz += sz[i]*sz[i];
  plhs[0]   = mxCreateDoubleMatrix(sumsz, 1, mxREAL);
  double *A = mxGetPr(plhs[0]);
  
  /* Invert all matrices. */
  getBlocksFromSparse(n, ia, ja, sa, A, nsz, sz);
}



int getBlocksFromSparse(int n, mwIndex *ia, mwIndex *ja, double *sa,  
                        double *A, int nsz, int *sz)
{
   
   int max_size = 0; 
   double *block  = A;
   
   int b;
   int row_begin = 0;
   
   for (b=0; b<nsz; ++b)
   {
      int i,j,k,r;
      int row_size = sz[b];
      int row_end = row_begin + row_size;

      /* Set to zero */
      for (i=0; i<row_size*row_size; ++i) block[i] = 0;
      
      for (r = row_begin; r<row_end; ++r)
      {
         for (k = ia[r]; k<ia[r+1]; ++k)
         {
           
            i                   = r     - row_begin;
            j                   = ja[k] - row_begin;
            block[i*row_size +j] = sa[k];
         }
      }
      row_begin = row_end;
      block = block + row_size*row_size;
   }
}
