#include <math.h>
#include <mex.h>
#include "fluid.h"

static int    nr;
static double *mobr;
static double *nw;
static double *no;
static double *srw;
static double *sro;

/*---------------------------------------------------------------------------*/
double fluxfun(double sw, const int reg)
/*---------------------------------------------------------------------------*/
{
  double den = 1 - srw[reg] - sro[reg];
  double so = (1 - sw - sro[reg])/den;
  sw = (sw - srw[reg])/den;
  sw = sw > 1.0 ? 1.0 : sw;
  sw = sw < 0.0 ? 0.0 : sw;

  double lw = pow(sw,nw[reg]);
  return lw/(lw + mobr[reg]*pow(so,no[reg]));
}

/*---------------------------------------------------------------------------*/
double dfluxfun(double sw, const int reg)
/*---------------------------------------------------------------------------*/
{
  double den = 1 - srw[reg] - sro[reg];
  double so = (1 - sw - sro[reg])/den;
  sw = (sw - srw[reg])/den;
  sw = sw > 1.0 ? 1.0 : sw;
  sw = sw < 0.0 ? 0.0 : sw;

  double lw    = pow(sw,nw[reg]);
  double lo    = pow(so,no[reg]);
  double denom = (lw + mobr[reg]*lo);

  return (mobr[reg]/den)*( nw[reg]*pow(sw,nw[reg]-1)*lo +
			   no[reg]*pow(so,no[reg]-1)*lw )/(denom*denom);
}

/*---------------------------------------------------------------------------*/
void init_fluid(const mxArray *arr)
/*---------------------------------------------------------------------------*/
{
  int i, n;

  /* Check sanity of input */
  if (! mxIsDouble(arr)) {
    mexErrMsgTxt("\nError in fluid object: parameter vector must be a vector of doubles\n");
    return;
  }
  n = (int) mxGetNumberOfElements(arr) - 2;
  if (n % 6 != 0) {
    mexErrMsgTxt("\nError in fluid object: parameter vector must contain 2+6*n elements\n");
    return;
  }

  /* Allocate sufficient memory to hold all regions */
  nr = n / 6;
  mobr = (double *) malloc( 5*nr * sizeof(double));
  nw   = mobr+nr;
  no   = nw+nr;
  srw  = no+nr;
  sro  = srw+nr;

  /* Extract data and initialize global variables */
  double *fprm = (double *) mxGetData(arr);
  for (i=0; i<nr; i++) {
    n = i*6;
    mobr[i] = (fprm[0]/fprm[1])*(fprm[7+n]/fprm[6+n]);
    nw[i]   = fprm[2+n];
    no[i]   = fprm[3+n];
    srw[i]  = fprm[4+n];
    sro[i]  = fprm[5+n];
    /* // mexPrintf("Fluid %d: %f %f %f %f %f\n", */
    /* //	      i+1, mobr[i], nw[i], no[i], srw[i], sro[i]); */
  }
  return;
}

/*---------------------------------------------------------------------------*/
void end_fluid()
/*---------------------------------------------------------------------------*/
{
  free(mobr);
  return;
}


/*---------------------------------------------------------------------------*/
int getNoSatRegions(void)
/*---------------------------------------------------------------------------*/
{
  return nr;
}


