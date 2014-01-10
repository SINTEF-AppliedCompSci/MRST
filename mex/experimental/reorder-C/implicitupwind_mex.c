/* Copyright 2011 (c) Jostein R. Natvig <Jostein.R.Natvig at sintef.no> */

#include <stdlib.h>
#include <string.h>
#include <mex.h>

#include "grid.h"
#include "fluid.h"
#include "mrst_api.h"
#include "twophasetransport.h"


/*
   To explicitly catch keyboard interrupt, a signal handler must be
   supplied, along with code to break off execution in case of user
   interrupt. Otherwise, Ctrl-C will not work in Matlab.
*/
int interrupt_signal = 0;
#include <signal.h>

static void
interruptHandler(int signal)
{
    interrupt_signal = signal;
}

void mexFunction(int nlhs,       mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{

    int i;
    struct UnstructuredGrid *grid;
    double *pv, *src, *darcyflux, *dt, *s0, *s;
    int     nsteps;
    int    *satnum;

    mwSize  m, n;

    /* return type for system call "signal" may not be typedef'ed on 
     * non-linux POSIX systems, such as osx */
    typedef void (*sighandler_t) (int);
    sighandler_t oldInterruptHandler;
    
   
    if (nrhs!=8) {
       mexErrMsgTxt("implicitupwind takes 8 input arguments");
    }
    if (nlhs>1) {
       mexErrMsgTxt("implicitupwind return a single argument");
    }
    grid = mrst_grid(prhs[0]);
    
    if (((int)mxGetNumberOfElements(prhs[1]) != grid->number_of_cells) ||
        ((int)mxGetNumberOfElements(prhs[2]) != grid->number_of_cells) ||
        ((int)mxGetNumberOfElements(prhs[3]) != grid->number_of_faces) ||
        ((int)mxGetNumberOfElements(prhs[5]) != grid->number_of_cells) ||
        ((int)mxGetNumberOfElements(prhs[6]) != grid->number_of_cells) )
    {
      mexErrMsgTxt("Dimensions mismatch in arguments...");
    }

    pv        = mxGetPr(prhs[1]);
    src       = mxGetPr(prhs[2]);
    darcyflux = mxGetPr(prhs[3]);
    dt        = mxGetPr(prhs[4]);
    nsteps    = mxGetNumberOfElements(prhs[4]);
    s0        = mxGetPr(prhs[5]);
    satnum    = mxGetData(prhs[6]);
    m         = mxGetM(prhs[5]);
    n         = mxGetN(prhs[5]);

    /* Supress compiler warning */
    (void) nrhs;

    init_fluid(prhs[7]);

    /* install handler for SIGINT signal */
    oldInterruptHandler = signal(SIGINT, interruptHandler);
    interrupt_signal = 0;

    plhs[0] = mxCreateDoubleMatrix(m, n, mxREAL);
    
    
    s = mxGetPr(plhs[0]);
    memcpy(s, s0, grid->number_of_cells * sizeof *s0);
    for (i=0; i<nsteps; ++i) {
       /* saturation at the beginning of the time step is replaced
        * by saturation at the end of the time step */
       twophasetransport(pv, src, dt[i], grid, darcyflux, satnum,
               s);
    }

    free_mrst_grid(grid);

    /* Restore old SIGINT handler */
    signal(SIGINT, oldInterruptHandler);

    if(interrupt_signal != 0)
    {
        mexErrMsgTxt("Interrupt signal caught in MRST mexFunction implicitupwind\n");
    }
}

/* Local Variables:    */
/* c-basic-offset:4    */
/* End:                */
