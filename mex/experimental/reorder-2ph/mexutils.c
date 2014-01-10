#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include <mex.h>
#include <string.h>
#include <stdio.h>
#include <stdarg.h>
#include "sparse.h"
#include "mexutils.h"



/******************************************************************************
 *                                                                            *
 * To explicitly catch keyboard interrupt, a signal handler must be           *
 * supplied, along with code to break off execution in case of user           *
 * interrupt. Otherwise, Ctrl-C will not work in Matlab.                      *
 *                                                                            *
 *****************************************************************************/


int interrupted  = 0;

#include <signal.h>



/*---------------------------------------------------------------------------*/
void interruptHandler(int signal_num)
/*---------------------------------------------------------------------------*/
{
  mexPrintf("\nsatsolver caught interrupt signal (%d)\n", signal_num);
  interrupted = 1;
}

/*---------------------------------------------------------------------------*/
void printSparse(sparse_t *S)
/*---------------------------------------------------------------------------*/
{
  /* Display the nonzero elements of the sparse array. */
  int m = S->m;
  int col, k;
  printf("%p\n", S);fflush(stdout);

  for (col=0; col<m; col++)
    {
      for (k = S->ia[col]; k < S->ia[col+1]; k++)
	mexPrintf("(%d,%d) = %g\n", S->ja[k], col, S->a[k]);
    }
}



/*
 * Copy Matlab sparse matrix to sparse_t matrix. Convert from mwIndex to index_t
 * NBNB!  Matlab uses compressed-sparse-columnm, while this code uses
 * compressed-sparse-row.
 */
/*---------------------------------------------------------------------------*/
sparse_t *getSparse(const mxArray *array)
/*---------------------------------------------------------------------------*/
{
  /* Get the starting positions of the data in the sparse array. */
  double  *a   = mxGetPr (array);
#if MATLABVERSION > 73
  mwIndex *ia  = mxGetJc (array);
  mwIndex *ja  = mxGetIr (array);
#else
  int     *ia  = mxGetJc (array);
  int     *ja  = mxGetIr (array);
#endif
  size_t   m   = mxGetN  (array);
  size_t   n   = mxGetM  (array);
  size_t   nnz = (size_t) ia[m];

  sparse_t *mat = sparse_alloc(m, n, nnz);

  size_t k;
  for (k=0; k<m+1; ++k) mat->ia[k] = (int) ia[k];
  for (k=0; k<nnz; ++k) mat->ja[k] = (int) ja[k];
  for (k=0; k<nnz; ++k) mat->a [k] = a[k];

  return mat;
}

/*---------------------------------------------------------------------------*/
mxArray* exportArray (const void *data, int ndata, int m, int n, mxClassID id)
/*---------------------------------------------------------------------------*/
{
  mxArray *arr;

  if (id == mxCHAR_CLASS)
  {
    const int max_size = 10000;
    if (strlen(data) > max_size)
    {
      mexErrMsgTxt("String input too large in exportArray. ");
    }
    arr = mxCreateString(data);
  }
  else
  {
    arr = mxCreateNumericMatrix(m, n, id, mxREAL);
    memcpy(mxGetPr(arr), data, ndata);
  }
  return arr;
}

/*---------------------------------------------------------------------------*/
void* copyMatlabVector(const mxArray *arr)
/*---------------------------------------------------------------------------*/
{
  void *v    = mxGetData (arr);
  int   n    = (int)mxGetNumberOfElements(arr);
  int   ne   = (int)mxGetElementSize(arr);

  void *copy = malloc (n * ne);

  memcpy (copy,  v, n * ne );
  return  copy;
}

/*---------------------------------------------------------------------------*/
void copyMatlabScalar(void *target, const mxArray *arr)
/*---------------------------------------------------------------------------*/
{
   void *v    = mxGetData    (arr);
   int   ne   = (int)mxGetElementSize(arr);
   int   n    = (int)mxGetNumberOfElements(arr);

   mxAssert (n==1, "Scalar argument expected");
   memcpy(target, v, ne);
}


/*---------------------------------------------------------------------------*/
void mxGetClassNameFromClassID(char *str, int n, mxClassID id)
/*---------------------------------------------------------------------------*/
{
  /* Fix this */

  switch(id)
  {
    case mxDOUBLE_CLASS:
      snprintf(str, n, "double");
      break;

    case mxINT32_CLASS:
      snprintf(str, n, "int32");
      break;

    default:
      snprintf(str, n, "some type");
  }
}



/*---------------------------------------------------------------------------*/
int copyMatlabData(void *target, int expected_id,
		   enum returnType type,  const mxArray *arr)
/*---------------------------------------------------------------------------*/
{
  void     *v  = mxGetData (arr);
  int       n  = (int)mxGetNumberOfElements(arr);
  mxClassID id = mxGetClassID(arr);
  int       ne = mxGetElementSize(arr);


  if (id != expected_id)
  {
    char str[64];
    mxGetClassNameFromClassID(str, 64, expected_id);
    mexPrintf   ("Found %s, expected %s.\n", mxGetClassName(arr), str);
    mexErrMsgTxt("Element type does not match expected type.");
  }

  if (type == VALUE)
  {
    memcpy(target, v, ne);
    return 1;
  }
  else /* (type == POINTER) */
  {
    void *copy;
    if (id == mxCHAR_CLASS)
    {
      /* Matlab string is vector of 2-byte elements. */
      copy = malloc ((n+1) * ne);
      mxGetString(arr, copy, n+1);
    }
    else
    {
      /* Assume that other datatypes have standard C type in data vector. */
      copy = malloc (n * ne);
      memcpy(copy, v, ne * n);
    }

    memcpy(target, &copy, sizeof(copy)); /* Copy pointer */
    return n;
  }
}




/*---------------------------------------------------------------------------*/
static int checkfields(int n, char *names[n], const mxArray *arr)
/*---------------------------------------------------------------------------*/
{
  int nfields    = mxGetNumberOfFields (arr);
  if (nfields != n)
    return 1;

  int i;
  for (i=0; i<n; ++i)
    if (NULL == mxGetField (arr,0,names[i]))
      return 1;

  return 0;
}

/*---------------------------------------------------------------------------*/
void  import_struct(void *out, int n, size_t num_elements[n], char *names[n],
		    int offset[n], mxClassID id[n],
		    enum returnType rtype[n], const mxArray *in)
/*---------------------------------------------------------------------------*/

{
  if(!mxIsStruct (in))
    mexErrMsgTxt ("Input to getiostruct must be a struct");

  int i;
  int nfields = mxGetNumberOfFields (in);

  if (checkfields(n, names, in))
    mexErrMsgTxt("Wrong number or wrong field names in struct.");

  for (i=0; i<nfields; ++i) {
    mxArray* field    = mxGetField (in, 0, names[i]);
    if (field == NULL) {
      char str[64];
      sprintf(str, "Field %s is empty!", names[i]);
      mexErrMsgTxt (str);
    }
    else
      num_elements[i] = copyMatlabData(out+offset[i], id[i], rtype[i], field);
  }
}
