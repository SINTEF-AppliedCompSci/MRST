/* Copyright 2011 (c) Jostein R. Natvig <Jostein.R.Natvig at sintef.no> */

#include <stdlib.h>

#ifdef MATLAB_MEX_FILE
#include <mex.h>
#include "grid.h"
#include "twophase.h"
#include "nlsolvers.h"
#include "fluid.h"
#else
#include <opm/core/grid.h>
#include <opm/core/transport/reorder/twophase.h>
#include <opm/core/transport/reorder/nlsolvers.h>
#include <opm/core/transport/reorder/fluid.h>
#endif

/* Parameters used in solution of single-cell boundary-value problem */
struct Parameters
{
    double s0;
    double influx;  /* sum_j min(v_ij, 0)*f(s_j) */
    double outflux; /* sum_j max(v_ij, 0)        */
    double dtpv;    /* dt/pv(i)                  */
    double div;
    int    region;
};


static struct Parameters get_parameters(struct SolverData *d, int cell);
static double residual(double s, void *data);
/* static double fluxfun(double s, int cell); */

void
destroy_solverdata(struct SolverData *d)
{
    if (d!=NULL)
    {
        free(d->fractionalflow);
    }
    free(d);
}

struct SolverData *
init_solverdata(struct UnstructuredGrid *grid, const double *darcyflux,
                const double *porevolume, const double *source,
                const int *satnum, double dt, double *saturation)
{
    int i;
    struct SolverData *d = malloc(sizeof *d);

    if(d!=NULL)
    {
        d->grid       = grid;
        d->darcyflux  = darcyflux;
        d->porevolume = porevolume;
        d->source     = source;
        d->satnum     = satnum;
        d->dt         = dt;

        d->saturation     = saturation;
        d->fractionalflow = malloc(grid->number_of_cells *
                                   sizeof *d->fractionalflow);
        if (d->fractionalflow == NULL)
        {
            destroy_solverdata(d);
            d = NULL;
        }
        for(i=0; i<grid->number_of_cells; ++i)
        {
            d->fractionalflow[i] = 0.0;
        }
    }
    return d;
}

/* Solver for single-cell bvp calls root-finder in nlsolvers.c */
void solvecell(void *data, struct NonlinearSolverCtrl *ctrl, int cell)
{
    struct SolverData   *d   = data;
    struct Parameters   prm = get_parameters(data, cell);
    
    d->saturation[cell] = find_zero(residual, &prm, ctrl);
    d->fractionalflow[cell] = fluxfun(d->saturation[cell], prm.region);
}


/* ====================== Internals =================================*/


/* static double  */
/* fluxfun(double s, int region) */
/* { */
/*     return s*s/(s*s + (1-s)*(1-s)); */
/* } */

/* Residual function r(s) for a single-cell bvp */
/*
 *     r(s) = s - s0 + dt/pv*(influx - outflux*f(s) )
 */
/* influx is water influx, outflux is total outflux */
static double
residual(double s, void *data)
{
    struct Parameters *p = data;
    return s - p->s0 +  p->dtpv*(-(p->div)*s+p->outflux*fluxfun(s, p->region) + p->influx);
}

static struct Parameters
get_parameters(struct SolverData *d, int cell)
{
    int i;
    struct UnstructuredGrid *g  = d->grid;
    struct Parameters        p;
    double flux;
    int f, other;

    p.div     = -(d->source[cell]);
    p.s0      = d->saturation[cell];
    p.influx  = d->source[cell] >  0 ? -d->source[cell] : 0.0;
    p.outflux = d->source[cell] <= 0 ? -d->source[cell] : 0.0;
    p.dtpv    = d->dt/d->porevolume[cell];
    if (d->satnum==NULL)
    {
        p.region =  0;
    }
    else
    {
        p.region = d->satnum[cell];
    }

    d->saturation[cell] = 0;
    for (i=g->cell_facepos[cell]; i<g->cell_facepos[cell+1]; ++i) {
        f = g->cell_faces[i];

        /* Compute cell flux*/
        if (cell == g->face_cells[2*f]) {
            flux  = d->darcyflux[f];
            other = g->face_cells[2*f+1];
        }
        else {
            flux  =-d->darcyflux[f];
            other = g->face_cells[2*f];
        }
	
	
        if (other != -1) {
            if (flux < 0.0) {
                p.influx  += flux*d->fractionalflow[other];
            }
            else {
                p.outflux += flux;
            }
	    p.div+= flux;
        }
    }
    return p;
}



/* Local Variables:    */
/* c-basic-offset:4    */
/* End:                */
