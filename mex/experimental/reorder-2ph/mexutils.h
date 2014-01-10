#ifndef _MEXUTILS_H_
#define _MEXUTILS_H_

#include <signal.h>

enum returnType {VALUE, POINTER};


void      interruptHandler (int signal_num);
void      printSparse      (sparse_t *S);
sparse_t *getSparse        (const mxArray *array);
void      import_struct    (void *out, int n, size_t num_elements[],
			    char *names[], int offset[],
			    mxClassID id[], enum returnType rtype[],
			    const mxArray *in);

mxArray  *exportArray      (const void *data, int ndata, int m, int n, mxClassID id);
void*     copyMatlabVector (const mxArray *arr);
void      copyMatlabScalar (void *target, const mxArray *arr);

#endif
