/*
 * Copyright 2010 (c) SINTEF ICT, Applied Mathematics.
 * Jostein R. Natvig <Jostein.R.Natvig at sintef.no>
 */



#include <mex.h>
#include <string.h>
#include "grid.h"
#include "../mrst_api/mrst_api.h"
#include "../mrst_api/call_umfpack.h"
#include "spu_implicit.h"



void
mexFunction(int nlhs,       mxArray *plhs[],
            int nrhs, const mxArray *prhs[])
{
    grid_t    *g;
    double    *s0, *s1, *s;
    double    *mob;
    double    *dmob;
    double    *dflux;
    double    *gflux;
    double    *src; 
    double    dt;
    double    norm = 0.0;
    
    double *tab;
    double h, x0;
    int *ntab;
    if (nrhs == 11 && nlhs == 2) {
        g       = mrst_grid(prhs[0]);
        s0      = mxGetPr  (prhs[1]);
        s1      = mxGetPr  (prhs[2]);
        h       = mxGetScalar(prhs[3]);
        x0      = mxGetScalar(prhs[4]);
        ntab    = (int*)mxGetData(prhs[5]);
        tab     = mxGetPr  (prhs[6]);
        dflux   = mxGetPr  (prhs[7]);
        gflux   = mxGetPr  (prhs[8]);
        src     = mxGetPr  (prhs[9]);
        dt      = mxGetScalar(prhs[10]);
        
        plhs[0] = mxCreateDoubleMatrix(g->number_of_cells, 1, mxREAL);
        s       = mxGetPr(plhs[0]);
        memcpy(s, s1, g->number_of_cells * sizeof *s);

        norm = spu_implicit(g, s0, s, h, x0, *ntab, tab, dflux, gflux, src, dt, callMWUMFPACK);

        plhs[1] = mxCreateDoubleScalar(norm);
        
        free_mrst_grid(g);
    }
    
    else {
        mexErrMsgTxt("Seven inputs, one output!\n");
    }
}
