/*===========================================================================
//
// File: mrst_fluid.c
//
// Authors: Ingeborg S. Ligaarden <Ingeborg.Ligaarden@sintef.no>
//          Jostein R. Natvig     <Jostein.R.Natvig@sintef.no>
//          Halvor M. Nilsen      <HalvorMoll.Nilsen@sintef.no>
//          Atgeirr F. Rasmussen  <atgeirr@sintef.no>
//          Bård Skaflestad       <Bard.Skaflestad@sintef.no>
//
//==========================================================================*/


/*
  Copyright 2011-2021 SINTEF Digital, Mathematics & Cybernetics.
*/

#include <stdlib.h>

#include <mex.h>

#include "mrst_fluid.h"

struct mrst_fluid_internal {
    mxArray *relperm;

    int      nphase;
    double  *mu;
    double  *rho;
};


/* ----------------------------------------------------------------------
 * Allocation, verification.
 * ---------------------------------------------------------------------- */
static int  verify_fluid          (const mxArray *FLUID);
static void fluid_internal_destroy(struct mrst_fluid_internal *i);

static struct mrst_fluid_internal *
fluid_internal_create(const mxArray *FLUID);

/* ----------------------------------------------------------------------
 * Declare fluid operations.
 * ---------------------------------------------------------------------- */
static int           nphase   (struct mrst_fluid *fluid);
static const double *viscosity(struct mrst_fluid *fluid);
static const double *density  (struct mrst_fluid *fluid);

static void          relperm  (struct mrst_fluid *fluid,
                               size_t             nsat ,
                               const double      *sat  ,
                               double            *kr   ,
                               double            *dkr  );

static void          destroy  (struct mrst_fluid *fluid);

/* ======================================================================
 * ====================================================================== */

/* ----------------------------------------------------------------------
 * Main public gateway to MRST fluid.
 * ---------------------------------------------------------------------- */
struct mrst_fluid *
mrst_fluid_create(const mxArray *FLUID)
/* ---------------------------------------------------------------------- */
{
    struct mrst_fluid *fluid;

    if (verify_fluid(FLUID)) {

        fluid = mxMalloc(1 * sizeof *fluid);

        if (fluid != NULL) {
            fluid->fluid_data = fluid_internal_create(FLUID);

            if (fluid->fluid_data != NULL) {
                fluid->nphase  = &nphase;
                fluid->visc    = &viscosity;
                fluid->dens    = &density;
                fluid->relperm = &relperm;
                fluid->destroy = &destroy;
            } else {
                destroy(fluid);
                fluid = NULL;
            }
        }
    } else {
        mexErrMsgTxt("Input FLUID must be a non-empty structure "
                     "containing fields\n\t'properties', and 'relperm'\n");
    }

    return fluid;
}


/* ----------------------------------------------------------------------
 * Implement memory allocation and fluid structure verification.
 * ---------------------------------------------------------------------- */

static int
verify_fluid(const mxArray *FLUID)
{
    int ok;

    ok =       mxIsStruct(FLUID) && !mxIsEmpty(FLUID);
    ok = ok && (mxGetFieldNumber(FLUID, "properties") >= 0);
    ok = ok && (mxGetFieldNumber(FLUID, "relperm"   ) >= 0);

    return ok;
}


static void
fluid_internal_destroy(struct mrst_fluid_internal *i)
{
    if (i != NULL) {
        if (i->relperm != NULL) { mxDestroyArray(i->relperm); }
        if (i->rho     != NULL) { mxFree        (i->rho    ); }
        if (i->mu      != NULL) { mxFree        (i->mu     ); }

        mxFree(i);
    }
}


static struct mrst_fluid_internal *
fluid_internal_create(const mxArray *FLUID)
{
    struct mrst_fluid_internal *i;
    const double *mu, *rho;

    int      nrhs, nlhs, status;
    size_t   p, np;
    mxArray *prhs[1];
    mxArray *plhs[2];

    nrhs    = 1;                /* fluid.properties */
    nlhs    = 2;                /* mu, rho */
    prhs[0] = mxGetField(FLUID, 0, "properties"); /* Get fhandle */

    /* [mu, rho] = feval(properties) % invoke fhandle */
    status = mexCallMATLAB(nlhs, plhs, nrhs, prhs, "feval");

    if (status == 0) {          /* OK */
        i = mxMalloc(1 * sizeof *i);

        np = mxGetN(plhs[0]);

        mxAssert (np > 0, "Number of phases must be positive.");
        mxAssert (mxGetM(plhs[0]) == 1,
                  "Viscosity must be a single row vector.");
        mxAssert (mxGetM(plhs[1]) == 1,
                  "Density must be a single row vector.");
        mxAssert (mxGetN(plhs[1]) == np,
                  "Density must be a row vector of the same "
                  "size as viscosity.");

        i->mu      = mxMalloc(np * sizeof *i->mu);
        i->rho     = mxMalloc(np * sizeof *i->rho);
        i->relperm = mxDuplicateArray(mxGetField(FLUID, 0, "relperm"));

        if ((i->mu == NULL) || (i->rho == NULL) || (i->relperm == NULL)) {
            fluid_internal_destroy(i);
            i = NULL;
        } else {
            mxAssert (mxIsDouble(plhs[0]), "Viscosity must be DOUBLE.");
            mxAssert (mxIsDouble(plhs[1]), "Density must be DOUBLE.");

            mu  = mxGetPr(plhs[0]);
            rho = mxGetPr(plhs[1]);

            for (p = 0; p < np; p++) {
                i->mu [p] = mu [p];
                i->rho[p] = rho[p];
            }

            i->nphase = (int) np;
        }

        /* Release resources acquired in properties() call */
        mxDestroyArray(plhs[1]);
        mxDestroyArray(plhs[0]);
    } else {
        i = NULL;
    }

    return i;
}


/* ----------------------------------------------------------------------
 * Fluid implementation.
 * ---------------------------------------------------------------------- */

static int
nphase(struct mrst_fluid *fluid)
{
    struct mrst_fluid_internal *i;

    mxAssert (fluid != NULL, "Huh: nphase() called without object?");

    i = fluid->fluid_data;

    return i->nphase;
}


static const double *
viscosity(struct mrst_fluid *fluid)
{
    struct mrst_fluid_internal *i;

    mxAssert (fluid != NULL, "Huh: visc() called without object?");

    i = fluid->fluid_data;

    return i->mu;
}


static const double *
density(struct mrst_fluid *fluid)
{
    struct mrst_fluid_internal *i;

    mxAssert (fluid != NULL, "Huh: dens() called without object?");

    i = fluid->fluid_data;

    return i->rho;
}


static void
relperm(struct mrst_fluid *fluid,
        size_t             nsat ,
        const double      *sat  ,
        double            *kr   ,
        double            *dkr  )
{
    struct mrst_fluid_internal *i;
    double                     *ptr;

    mwSize   s, ns, p, np;
    int      nlhs   ,  nrhs   , status;
    mxArray *plhs[2], *prhs[2];

    mxAssert (fluid != NULL, "Huh: relperm() called without object?");

    i  = fluid->fluid_data;
    ns = nsat;
    np = i->nphase;

    nrhs    = 2;
    prhs[0] = i->relperm;
    prhs[1] = mxCreateDoubleMatrix(ns, np, mxREAL);

    /* Transpose input for natural access in M */
    for (p = 0, ptr = mxGetPr(prhs[1]); p < np; p++) {
        for (s = 0; s < ns; s++) {
            *ptr++ = sat[s*np + p];
        }
    }

    nlhs = 1 + (dkr != NULL);

    /*  kr       = feval(relperm, s) % or
     * [kr, dkr] = feval(relperm, s) */
    status = mexCallMATLAB(nlhs, plhs, nrhs, prhs, "feval");

    if (status == 0) {          /* OK */
        /* Extract output.  Transpose for natural access in C */
        mxAssert (mxGetM(plhs[0]) == ns,
                  "Rel-perm must be computed at each saturation point.");
        mxAssert (mxGetN(plhs[0]) == np,
                  "Rel-perm must be computed for each phase.");

        for (p = 0, ptr = mxGetPr(plhs[0]); p < np; p++) {
            for (s = 0; s < ns; s++) {
                kr[s*np + p] = *ptr++;
            }
        }

        if (dkr != NULL) {
            mxAssert (mxGetM(plhs[0]) == ns,
                      "Rel-perm derivative must be computed at "
                      "each saturation point.");
            mxAssert (mxGetN(plhs[1]) == np * np,
                      "Rel-perm derivative must be np-by-np.");

            for (p = 0, ptr = mxGetPr(plhs[1]); p < np * np; p++) {
                for (s = 0; s < ns; s++) {
                    dkr[s*np*np + p] = *ptr++;
                }
            }

            mxDestroyArray(plhs[1]); /* Release 'dkr' */
        }

        mxDestroyArray(plhs[0]);     /* Release 'kr' */
    }

    mxDestroyArray(prhs[1]);         /* Release 's' */
}


static void
destroy(struct mrst_fluid *fluid)
{
    if (fluid != NULL) {
        fluid_internal_destroy(fluid->fluid_data);

        mxFree(fluid);
    }
}
