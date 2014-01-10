/*
 * Copyright 2010 (c) SINTEF ICT, Applied Mathematics.
 * Jostein R. Natvig <Jostein.R.Natvig at sintef.no>
 */



#include <mex.h>
#include "grid.h"
#include "mrst_api.h"
#include "spu_explicit.h"

/* ---------------------------------------------------------------------- */
void
mexFunction(int nlhs,       mxArray *plhs[],
            int nrhs, const mxArray *prhs[])
/* ---------------------------------------------------------------------- */
{
    grid_t    *g;
    double    *s0, *s;
    double    *mob;
    double    *dflux;
    double    *gflux;
    double    *src; 
    double    dt;

 
    if (nrhs == 7 && nlhs == 1) {
        g       = mrst_grid(prhs[0]);
        s0      = mxGetPr  (prhs[1]);
        mob     = mxGetPr  (prhs[2]);
        dflux   = mxGetPr  (prhs[3]);
        gflux   = mxGetPr  (prhs[4]);
        src     = mxGetPr  (prhs[5]);
        dt      = mxGetScalar(prhs[6]);
        
        plhs[0] = mxCreateDoubleMatrix(g->number_of_cells, 1, mxREAL);
        s       = mxGetPr(plhs[0]);
        
        spu_explicit(g, s0, s, mob, dflux, gflux, src, dt);
        
        free_mrst_grid(g);
    }
    
    else {
        mexErrMsgTxt("Seven inputs, one output!\n");
    }
}
