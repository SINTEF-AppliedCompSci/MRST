#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <math.h>

#include <mex.h>

#include "sparse.h"
#include "nlsolvers.h"
#include "system.h"
#include "fluid.h"


static void compute_upwind_flux(sparse_t*, double*, double*, sparse_t**, double *work);
static void extractScalarSubProblem(int, system_t*, double*);
static enum Method {RIDDER, REGULAFALSI, BISECTION} method = REGULAFALSI;


static double G(double s, void *data)
/*---------------------------------------------------------------------------
  ---------------------------------------------------------------------------*/
{
  double *D = (double *) data;
  return s + D[0]*fluxfun(s, (int) D[3]) - D[1];
}

static double dG (double s, void *data)
/*---------------------------------------------------------------------------
  ---------------------------------------------------------------------------*/
{
  double *D = (double *) data;
  return 1.0 + D[0]*dfluxfun(s,(int) D[3]);
}


int solveScalar(int cell, void *d, void *r)
/*---------------------------------------------------------------------------
 cell  - cell number
 d     - data structure representing system (flux, pore volume, old saturations, ...)
 r     - report  (number of iterations, ...)
 opt   - options (tolerance, max iterations, solver type, ...)
 ---------------------------------------------------------------------------*/
{
  system_t *sys = d;
  report_t *rep = r;

  double       TOL    = sys->tolerance;/* 1.0e-6; */
  const int    maxit  = sys->maxit;
  double      *s      = sys->u;

  double       data[4];

  extractScalarSubProblem(cell, sys, data);
  data[2] = sys->u_old[cell];                       /* // initial guess */
  data[3] = (sys->satnum) ? sys->satnum[cell] : 0;  /* // saturation region */

  /* Check if solution in [0,1] exist */
  if (( data[0] + 1 < data[1] ) ||  ( data[1] < 0 )) {
    if  ( data[1] <= 0.0)
      s[cell] = 0.0;
    else
      s[cell] = 1.0;
  }
  else {
    switch (method) {
      default:
      case BISECTION:
	s[cell] = bisection  ( &G, data, TOL, maxit, &rep->iterations[cell]); break;
      case RIDDER:
	s[cell] = ridder     ( &G, data, TOL, maxit, &rep->iterations[cell]); break;
      case REGULAFALSI:
	s[cell] = regulafalsi( &G, data, TOL, maxit, &rep->iterations[cell]); break;
    }
  }

  sys->f[cell]   = fluxfun(s[cell], (int) data[3]);
  return 0;
}


static void extractScalarSubProblem(int cell, system_t *sys, double *D)
/*---------------------------------------------------------------------------
  ---------------------------------------------------------------------------*/
{
  int k;
  double sumIn  = 0.0;
  double sumOut = 0.0;

  for (k=sys->V->ia[cell]; k<sys->V->ia[cell+1]; ++k)
  {
    int j = sys->V->ja[k];
    assert ( j <= cell );

    if (j<cell) /* // Incoming flux */
      sumIn -= sys->V->a[k]*sys->f[j];
    else        /* // Outgoing flux and negative source *\/ */
      sumOut   += sys->V->a[k];
  }

  /* Incoming source */
  sumIn -= sys->Qp[cell]*fluxfun(1.0, sys->satnum[cell]);

  D[0] = sumOut *sys->dt;
  D[1] = sumIn  *sys->dt + sys->u_old[cell];

  return;
}


void init_system(system_t *sys, double dt, sparse_t *flux, double *sources,
		 double *volumes, double *u, int *satnum,
		 char *solver, int maxit, double tolerance, double   *work)
/*---------------------------------------------------------------------------
  work must be 2n * sizeof(double) bytes long.  Work is here to avoid
  malloc/free inside this function.
  ---------------------------------------------------------------------------*/
{
  sparse_t *upwind_flux = {0};

  if (strstr(solver,"ridder"))
    method = RIDDER;
  else if (strstr(solver,"regulafalsi"))
    method = REGULAFALSI;
  else if (strstr(solver,"bisection"))
    method = BISECTION;
  else {
    mexPrintf("Error: nonexisting nonlinear solver (%s)\n", solver);
    return;
  }

  /* Allocates a new sparse matrix. */
  compute_upwind_flux(flux, sources, volumes, &upwind_flux, work);

  sys->n      = (size_t) flux->m;
  sys->dt     = dt;
  sys->V      = upwind_flux;

  sys->u      = u;
  sys->u_old  = work;
  sys->f      = work +   sys->n;
  sys->Qp     = sources;

  if (satnum != NULL) {
    int i;
    int nr = getNoSatRegions()-1;
    for (i=0; i<sys->n; i++)
      if ((satnum[i]<0) || (satnum[i] > nr)) {
	mexErrMsgTxt("\nError in init_system: invalid saturation region\n");
	return;
      }
  }
  sys->satnum = satnum;

  memcpy(sys->u_old, sys->u, sys->n*sizeof(double));

  sys->maxit = maxit;
  sys->tolerance = tolerance;

  return;
}

static void compute_upwind_flux(sparse_t *Fin, double *source, double *volume,
				sparse_t **upwind_flux, double *work)
/*---------------------------------------------------------------------------
  ---------------------------------------------------------------------------*/
{
  /* work must be n*sizeof(double) bytes long. */
  int       i,k;
  const int m = Fin->m;
  const int N = Fin->ia[m];

  sparse_t *F = sparse_alloc(m, m, N+m);

  /* Copy Fin to F and divide by PV */
  for (i=0; i<m+1; ++i)
    F->ia[i] = Fin->ia[i]+i;
  for (i=0; i<m; ++i)
    for (k=Fin->ia[i]; k<Fin->ia[i+1]; ++k) {
      F->ja[k+i] =  Fin->ja[k];
      F->a [k+i] = -Fin->a [k]/volume[i];
    }

  for (i=0; i<N; ++i) work[i] = 0.0;

  /* Add outgoing fluxes */
  double *r = Fin->a;
  int    *q = Fin->ja;
  for (i=0; i<N; ++i)
    work[*q++] += *r++;

  /* Add negative source */
  for (i=0; i<m; ++i) {
    work  [i] += (source[i] < 0) ? -source[i]               : 0.0;
    source[i]  = (source[i] > 0) ? -source[i]*1.0/volume[i] : 0.0;
  }

  /* Put outgoing fluxes on diagonal. */
  for (i=0; i<m; ++i) {
    k = F->ia[i+1]-1;
    F->ja[k] = i;
    F->a [k] = work[i]/volume[i];
  }

  *upwind_flux = F;
}

