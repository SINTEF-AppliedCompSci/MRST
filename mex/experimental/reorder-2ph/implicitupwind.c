/*------------------------------------------------------------
 File       : implicitupwind.c
 Created    : Tue, 15 Nov  2011
 Author     : Jostein R. Natvig <Jostein.R.Natvig@sintef.no>
            : Knut-Andreas Lie  <Knut-Andreas.Lie@sintef.no>

 Description
 Mex interface to reordered first-order upwind solver without
 compressibility. The solver assumes that there are no loops
 in the velocity field.

------------------------------------------------------------*/
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#define ENABLE_MEXERRMSGTXT 1
#define CHECK_DIVERGENCE    0

#include <string.h>
#include <mex.h>
#include <matrix.h>
#include <math.h>
#include <stdio.h>
#include <stdarg.h>

#include "fluid.h"
#include "sparse.h"
#include "tarjan.h"
#include "system.h"
#include "utils.h"
#include "mexutils.h"

/*#include "today.h"*/

extern int    interrupted;


/* ---- Settings ----------------------------- */
struct options
{
  int    max_iterations;
  int    substeps;
  double tolerance;
  int    verbosity;
  char   *scalar_solver;
};


/* ----- Statistics and diagnostic output ------*/
static const int strlength = 16384;
struct report
{
  double average_iterations;
  int    max_iterations;
  int    num_zeroiterations;
  double internal_time;
  char   message[16384];
};

/* Declaration of internal functions */
static int  isFull      (const mxArray *);
static int  checkInput  (const int, const mxArray**);
static void get_options (struct options *, size_t, size_t[], const mxArray *);
static int  solve       (double, sparse_t*, double*, double*, double*, int *,
			 struct options*, struct report*);
static void export_report(const struct report *, mxArray **);
static int  Message (char **, int *, const char *, ...);

#if CHECK_DIVERGENCE
static int  checkDivergence (const sparse_t*, const double*);
#endif


void mexFunction(int nlhs,       mxArray *plhs[],
		 int nrhs, const mxArray *prhs[])
/*---------------------------------------------------------------------------
   The mex interface to the solver
  ---------------------------------------------------------------------------*/
{

  if (checkInput(nrhs, prhs)) return;

  /* Install handler for SIGINT signal */
  sighandler_t oldInterruptHandler = signal(SIGINT, interruptHandler);
  interrupted    = 0;

  /* Initialize fluid object */
  init_fluid(prhs[5]);

  /* Fetch input and allocate output for solver */
  double   tf;         copyMatlabScalar(&tf,prhs[3]);
  sparse_t *flux     = getSparse        (prhs[0]);
  double   *pv       = copyMatlabVector (prhs[1]);
  plhs[0]            = mxDuplicateArray (prhs[2]);
  double   *sat      = mxGetPr          (plhs[0]);
  double   *q        = copyMatlabVector (prhs[4]);
  int      *satnum   = (mxIsEmpty(prhs[6])) ? NULL : copyMatlabVector (prhs[6]);

  /* Read options that control simulation */
  struct options opt = {0};
  size_t fieldsize[5];
  get_options(&opt, sizeof fieldsize / sizeof fieldsize[0], fieldsize, prhs[7]);

  /* Call solver */
  struct report r;
  solve(tf, flux, pv, q, sat, satnum, &opt, &r);

  /* Free local copies of Matlab data */
  free(q);
  free(pv);
  sparse_free(flux);
  free(satnum);

  end_fluid();

  if (nlhs == 2)
    export_report (&r, &plhs[1]);

  /* Restore old SIGINT handler */
  signal(SIGINT, oldInterruptHandler);

#if ENABLE_MEXERRMSGTXT
  /* Force an interrupt in Matlab */
  if ( interrupted )
    mexErrMsgTxt("User interrupt");
#endif
}


/*==========================================================
 *             INTERNAL FUNCTIONS
 ==========================================================*/

static int solve(double tf, sparse_t *flux, double *volume, double *source,
		 double *u, int *satnum, struct options *opt, struct report  *r)
/*---------------------------------------------------------------------------
 * Input arguments:
 *
 *     tf       - final time
 *     dt       - magnitude of time step
 *     flux     - nxn sparse matrix of (positive) in-fluxes.
 *                (Confusing: The same as Fin in computeFlux)
 *     volume   - n-vector holding pore volumes
 *     sources  - n-vector holding sources and sinks for each cell
 *     u        - n-vector holding current satruations S.
 *     satnum   - n-vector holding index to saturation regions
 *
 * Output:
 *     u        - n-vector of new saturations S.
 *
 * The implicit upwind scheme for incompressible two-phase transport
 *
 *    s^n - s^(n-1) + dt/volume*flux*f(s) = sources
 *
 * is used to advance the solution from t=0 to t=tf in time steps of
 * size dt.
 ---------------------------------------------------------------------------*/
{
  int    i, begin, end;
  double t;

  char *str = r->message;
  int  len  = strlength;

  const int nc = flux->m;
  int      *rowP    = malloc (  nc    * sizeof(int));
  int      *P       = malloc ( (nc+1) * sizeof(int));
  double   *work    = malloc (  5*nc  * sizeof(double));

  /* Reorder unknowns and permute system */
  i = tarjan(flux->m, flux->ia, flux->ja, rowP, P, (int *) work);
  tic(0);
  if (nc != i) {
    mexPrintf("Solver not constructed to handle loops\n"
	      "Number of components (%d) different than number of cells (%d)", i, nc);
    return 1;
  }
  sparse_permute(flux, rowP);
  memcpy(work, volume, nc*sizeof(double)); for(i=0; i<nc; ++i) volume[i] = work[rowP[i]];
  memcpy(work, source, nc*sizeof(double)); for(i=0; i<nc; ++i) source[i] = work[rowP[i]];
  memcpy(work, u,      nc*sizeof(double)); for(i=0; i<nc; ++i) u     [i] = work[rowP[i]];
  int *buf = malloc( nc*sizeof(int));
  memcpy(buf, satnum, nc*sizeof(int)); for(i=0; i<nc; ++i) satnum[i] = buf[rowP[i]];
  free(buf);

  /* Write report */
  if(opt->verbosity>0){
    if (!sparse_is_triangular(flux, P))
      mexWarnMsgTxt ("The permuted matrix of fluxes is not triangular\n");

    Message (&str, &len, "\n---------------------------------------------------\n");
    Message (&str, &len, "Reordering            :    %2.2e sec\n", toc(0));
    Message (&str, &len, "Nonlinear solver      :    %s\n",    opt->scalar_solver);
    Message (&str, &len, "Saturation regions    :    %d\n",   getNoSatRegions());
    Message (&str, &len, "Number of fluxes      :    %d\n",   (int)flux->ia[flux->m]);
    Message (&str, &len, "Number of cells       :    %d\n",   flux->m);
    Message (&str, &len, "Subdomains            :    %d\n",   nc);
    Message (&str, &len, "Max iterations        :    %d\n",   opt->max_iterations );
    Message (&str, &len, "Tolerance             :    %2.2e\n", opt->tolerance);
    Message (&str, &len, "Computing solution    :    ");  fflush(stdout);
  }

  /* Initialize system */
  system_t sys = {0};
  init_system(&sys, tf/opt->substeps, flux, source, volume, u, satnum,
	      opt->scalar_solver, opt->max_iterations, opt->tolerance, work);

  /* Init report */
  report_t rep = {0};
  rep.iterations = calloc(flux->m, sizeof(int));

  /* Main loop */
  t  = 0.0;                /* physical time */
  r->internal_time = 0.0;  /* cpu time */
  while ((t < tf) && !interrupted)
  {
    if(opt->verbosity>0) {
      if (t>0) Message(&str, &len, "                      :    ");
      Message (&str, &len,"(dt=%2.2e, runtime:", tf/opt->substeps);
    }
    tic(0);

    /* Compute solution */
    for (i=0; i<nc; ++i) {
      begin  = P[i];
      end    = P[i+1];
      if (end-begin > 1) {
	mexPrintf("Solver not constructed to handle loops\n"
		 "Loop in component %d (%d...%d)\n", i, begin, end-1);
	return 1;
      }

      if (solveScalar(begin, &sys, &rep)) {
	mexPrintf("Error occured in component %d (%d...%d)\n", i, begin, end-1);
	return 1;
      }
      if (interrupted) { mexPrintf("interrupted"); return 1; }
    }

    /* Swap u and u_old */
    double *tmp = sys.u; sys.u = sys.u_old; sys.u_old = tmp;

    /* Advance time */
    t  = t + sys.dt;
    sys.dt = (sys.dt < tf - t) ? sys.dt : tf - t;

    /* Add to timer */
    double elapsed = toc(0);
    if(opt->verbosity>0) Message(&str, &len, "%2.2e sec)\n", elapsed);
    r->internal_time += elapsed;
  }

  if (u == sys.u)
    memcpy(u, sys.u_old, sys.n*sizeof(double));

  /* Permute final solution back to original sequence */
  memcpy(work, u, nc*sizeof(double));
  for(i=0; i<nc; ++i) u[rowP[i]] = work[i];

  int max=0;
  int sum=0;
  int num0=0;
  for (i=0; i<nc; ++i) {
    sum += rep.iterations[i];
    max  = (max < rep.iterations[i]) ? rep.iterations[i] : max;
    if (rep.iterations[i]==0) num0++;
  }
  if(opt->verbosity>0){
    Message (&str, &len, "Average iterations    :    %2.2f\n", sum*1.0/nc);
    Message (&str, &len, "Max iterations        :    %d\n", max);
    Message (&str, &len, "Number of zero it.    :    %d\n", num0);
    Message (&str, &len, "---------------------------------------------------\n"); fflush(stdout);
  }
  r->average_iterations = sum*1.0/nc;
  r->max_iterations     = max;
  r->num_zeroiterations = num0;

  /* Free report */
  free(rep.iterations);

  sparse_free(sys.V);
  free (rowP);
  free (P);
  free (work);

  return 0;
}



static void export_report(const struct report *in, mxArray **out)
/*---------------------------------------------------------------------------
  Write report about number of iterations etc
  ---------------------------------------------------------------------------*/
{
  const mwSize  nfields = 5;
  const char   *fieldnames[5] = {
    "average_iterations", "max_iterations", "num_zeroiterations",
    "internal_time", "message"};

  *out = mxCreateStructMatrix(1, 1, nfields, fieldnames);
  mxSetField (*out, 0, "average_iterations",
	      exportArray(&in->average_iterations, sizeof(double), 1, 1, mxDOUBLE_CLASS));
  mxSetField (*out, 0, "max_iterations",
	      exportArray(&in->max_iterations,     sizeof(int),    1, 1, mxINT32_CLASS));
  mxSetField (*out, 0, "num_zeroiterations",
	      exportArray(&in->num_zeroiterations, sizeof(int),    1, 1, mxINT32_CLASS));
  mxSetField (*out, 0, "internal_time",
	      exportArray(&in->internal_time,      sizeof(double), 1, 1, mxDOUBLE_CLASS));
  mxSetField (*out, 0, "message",
	      exportArray(in->message, sizeof(char), 1, strlen(in->message), mxCHAR_CLASS));

  return;
}



static void get_options (struct options *out, size_t nelms,
			 size_t *num_elements, const mxArray *in)
/*---------------------------------------------------------------------------
  Extract the variable options that are used to govern the solver
    - maximum number of iterations
    - number of substeps used in the time stepping
    - verbosity of output
    - name of the nonlinear scalar solver
  ---------------------------------------------------------------------------*/
{
  const int n = 5;
  int  offset[5];
  char *names [5] = {"max_iterations", "substeps", "tolerance",
		       "verbosity","scalar_solver"};
  enum returnType rtype[5] = {VALUE, VALUE, VALUE, VALUE,POINTER};

  offset[0] = offsetof(struct options, max_iterations);
  offset[1] = offsetof(struct options, substeps);
  offset[2] = offsetof(struct options, tolerance);
  offset[3] = offsetof(struct options, verbosity);
  offset[4] = offsetof(struct options, scalar_solver);

  mxClassID  id[5]= {mxINT32_CLASS, mxINT32_CLASS, mxDOUBLE_CLASS,
		     mxINT32_CLASS, mxCHAR_CLASS};
  mxAssert ((size_t) n == nelms,"Disparate number of fields.");
  import_struct(out, n, num_elements, names, offset, id, rtype, in);
  return;
}



static int Message (char **buf, int *bufsize, const char *format, ...)
/*---------------------------------------------------------------------------
  Displaying message to screen and storing in a buffer
  ---------------------------------------------------------------------------*/
{
  int size;

  va_list ap;
  va_start(ap, format);

  size = vsprintf(*buf, format, ap);
  mexPrintf("%s", *buf);
  *buf     += size;
  *bufsize += size;
  va_end(ap);

  return size;
}


static int checkInput(const int nrhs, const mxArray *prhs[])
/*----------------------------------------------------------------------
  Check that we are getting sensible input from Matlab
  ---------------------------------------------------------------------*/
{
  if (nrhs != 8) {
    mexErrMsgTxt("Error: please call me with EIGHT input arguments\n");
    return 1;
  }

  int m = mxGetM(prhs[0]);
  int n = mxGetN(prhs[0]);

  /* Flux matrix */
  if (!mxIsSparse(prhs[0]) || (m!=n))  {
    mexErrMsgTxt("First argument must be a square sparse matrix.\n");
    return 1;
  }

  /* Pore volume */
  if ((int)mxGetNumberOfElements(prhs[1]) != m || !isFull(prhs[1])) {
    mexErrMsgTxt(
	  "\nInput error:\n"
	  "Argument 2 (porevolume) must have the same number of elements as\n"
	  "the number of rows in argument 1.\n" );
    return 1;
  }

  /* Initial saturation */
  if ((int)mxGetNumberOfElements(prhs[2]) != m || !isFull(prhs[2])) {
    mexErrMsgTxt(
	  "\nInput error:\n"
	  "Argument 3 (initial saturation) must have the same number of elements\n"
	  "as the number of rows in argument 1.\n");
    return 1;
  }

  /* Final time */
  if (mxGetNumberOfElements(prhs[3]) != 1 || !mxIsDouble(prhs[3])){
    mexErrMsgTxt(
	  "\nInput error:\n"
	  "Argument 4 (final time) must be a double.\n" );
    return 1;
  }

  /* Source term */
  if ((int)mxGetNumberOfElements(prhs[4]) != m || !isFull(prhs[4])) {
    mexErrMsgTxt(
	  "\nInput error:\n"
	  "Argument 5 (source term) must have the same number of elements as\n"
	  "the number of rows in argument 1.\n");
    return 1;
  }

  /* Fluid parameters */
  if (! mxIsDouble(prhs[5])) {
    mexErrMsgTxt(
	   "\nInput error:\n"
	   "Argument 6 (fluid parameters) must be a vector of doubles");
    return 1;
  }

  /* Saturation regions */
  if (!( mxIsEmpty(prhs[6]) ||
	 ( mxIsInt32(prhs[6]) && ((int)mxGetNumberOfElements(prhs[6])==m)) )) {
    mexErrMsgTxt(
           "\nInput error:\n"
	   "Argument 7 (satnum) must either be [] or an int32 vector with\n"
	   "the same number of elements as the number of rows in argument 1.\n");
    return 1;
  }

  /* Options */
  if ((int)mxGetNumberOfElements(prhs[7]) != 1 || !mxIsStruct(prhs[7])) {
    mexErrMsgTxt(
	  "\nInput error:\n"
	  "Argument 8 (options) must have one element and be a struct.\n" );
    return 1;
  }

  return 0;
}



int isFull(const mxArray *a)
/*---------------------------------------------------------------------------
  Check if we have a full matrix
  ---------------------------------------------------------------------------*/
{
  return !(mxIsSparse(a) || mxIsStruct(a) || mxIsCell(a));
}



#if CHECK_DIVERGENCE
static int checkDivergence(const sparse_t *A, const double *Q)
/*---------------------------------------------------------------------------
  Check if the flux field is divergence free (as is assumed in the solver)
  ---------------------------------------------------------------------------*/
{
  int     m = A->m;
  int     i,j,k;
  double  val;

  double  *rowsum = calloc(m, sizeof(double));

  for(i=0; i<m; ++i)
  {
    for (k=A->ia [i]; k<A->ia [i+1]; ++k)
    {
      j   = A->ja [k];
      val = A->a  [k];

      if (val < 0.0)
      {
	mexErrMsgTxt ("Error!\n\n\tPositive fluxes expected.\n\n");

	free (rowsum);
	return 1;
      }
      rowsum [j] -= val;
      rowsum [i] += val;
    }
    rowsum [i] += Q[i];
  }

  double maxflux = 0.0;
  for (k=0; k<A->ia[m]; ++k)
    if (fabs (A->a [k]) > maxflux)
      maxflux = fabs (A->a [k]);

  double infnorm = 0.0;
  for(i=0; i<m; ++i)
    if (fabs (rowsum [i]) > infnorm)
      infnorm = fabs (rowsum [i]);

  free(rowsum);
  if (infnorm/maxflux > 1e-6) {
    mexErrMsgTxt("Error!\n\n\t||div(v)||/||v|| > 1e-6 (%e)\n\n", infnorm);
    return 1;
  }
  else
    message(2, "OK\n");

  return 0;
}
#endif
