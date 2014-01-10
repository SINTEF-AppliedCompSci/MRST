/* Copyright 2011 (c) Jostein R. Natvig <Jostein.R.Natvig at sintef.no> */

#include <stdlib.h>
#include <string.h>
#include <mex.h>
#include <matrix.h>

#include "grid.h"
#include "fluid.h"
#include "mrst_api.h"
#include "reordersequence.h"


/*
   To explicitly catch keyboard interrupt, a signal handler must be
   supplied, along with code to break off execution in case of user
   interrupt. Otherwise, Ctrl-C will not work in Matlab.
*/
int interrupt_signal  = 0;
#include <signal.h>

static void
interruptHandler(int signal)
{
    interrupt_signal = signal;
}

void mexFunction(int nlhs,       mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{

    struct UnstructuredGrid *grid;

    double  *darcyflux;
    int      i;
    mwSize   sz;
    mxArray *mxsequence;
    mxArray *mxcomponents;
    int     *sequence;
    int     *components;
    int      ncomponents;

    mwSize   dims[2];
    mwSize   ndims = 2;

    /* return type for system call "signal" may not be typedef'ed on
     * non-linux POSIX systems, such as osx */
    typedef void (*sighandler_t) (int);
    sighandler_t oldInterruptHandler;
    
    if (nrhs!=2)
    {
       mexErrMsgTxt("reordersequence takes two input arguments");
    }
    grid = mrst_grid_topology(prhs[0]);
    sz   = grid->number_of_cells;
    
    if (grid==NULL)
    {
       mexErrMsgTxt("Invalid grid as first argument");
    }
    if (grid->number_of_faces != (int)mxGetNumberOfElements(prhs[1]))
    {
       mexErrMsgTxt("Number of fluxes (in second argument) must equal number of faces");
    }
    darcyflux = mxGetPr(prhs[1]);
    
    /* install handler for SIGINT signal */
    oldInterruptHandler = signal(SIGINT, interruptHandler);
    interrupt_signal    = 0;

    
    /* create arrays in which to return sequence */
    mxsequence   = mxCreateNumericMatrix(sz,   1, mxINT32_CLASS, mxREAL);
    mxcomponents = mxCreateNumericMatrix(sz+1, 1, mxINT32_CLASS, mxREAL);
    
    sequence     = mxGetData(mxsequence);
    components   = mxGetData(mxcomponents);
    
    
    /* compute sequence using Tarjans algorithm */
    compute_sequence(grid, darcyflux, sequence, components, &ncomponents);
    
    
    /* in Matlab, indexing starts with one... */
    for (i=0; i<grid->number_of_cells; ++i)
    {
       sequence[i] += 1;
    }

    
    /* return sequence and optionally, components */
    plhs[0] = mxsequence;
    if (nlhs == 1)
    {
       mxDestroyArray(mxcomponents);
    }
    
    if(nlhs==2)
    {
       if (ncomponents != (int)sz)
       {
          components = mxRealloc(components, ncomponents);
          mxSetData(mxcomponents, components);
          dims[0]   = ncomponents;
          dims[1]   = 1;
          mxSetDimensions(mxcomponents, dims, ndims);
       }

       /* in Matlab, indexing starts with one... */
       for (i=0; i<ncomponents+1; ++i)
       {
          components[i] +=1;
       }
       plhs[1] = mxcomponents;
    }

    
    /* clean up */
    free_mrst_grid(grid);


    /* Restore old SIGINT handler */
    signal(SIGINT, oldInterruptHandler);

    if(interrupt_signal != 0)
    {
        mexErrMsgTxt("Interrupt signal caught in MRST mexFunction\n");
    }
}

/* Local Variables:    */
/* c-basic-offset:4    */
/* End:                */
