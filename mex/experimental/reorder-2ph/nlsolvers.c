/*------------------------------------------------------------
 File       : nlsolvers.c
 Created    : Thu Oct 16 09:38:02 2008
 Author     : Jostein R. Natvig <Jostein.R.Natvig@sintef.no>
 Revision   : $Id$


 Description


 Changes

------------------------------------------------------------*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include "nlsolvers.h"


#if MATLAB_MEX_FILE
#include <mex.h>
extern int interrupted;
#define print   mexPrintf
#define malloc  mxMalloc
#define calloc  mxCalloc
#define realloc mxRealloc
#define free    mxFree
#endif


static const char no_root_str[]=
  "  In %s:\n"
  "  With G(0) =% 5f, G(1) =% 5f, G(s) cannot have a zero in [0,1]!\n";

static const double EPS = 1e-14;


/* Start with bracket [0,1] with G(0)*G(1)<0.  */
/*---------------------------------------------------------------------------*/
double
ridder (double (*G)(double, void*), void *data, double tol, int maxit, int *iterations)
/*---------------------------------------------------------------------------*/
{
  *iterations = 0;

  double s2 = ((double *) data)[2];
  double G2 = G(s2, data);
  if (fabs(G2)<tol) return s2;

  /* Initial guess is interval [0,1] of course */
  double s0=0.0;
  double G0=G(s0, data);
  if (fabs(G0)<EPS) return s0;

  double s1=1.0;
  double G1=G(s1, data);
  if (fabs(G1)<EPS) return s1;

  if (G0*G1 > 0)
  {
    print(no_root_str, "ridder", G0, G1);
    return -1.0;
  }

  if (G0>G1)
  {
    double swap;
    swap = G0;
    G0   = G1;
    G1   = swap;
  }

  double s3=0;
  double G3=10;

  int it = 0;
  while ( (fabs(G3) > tol)  &&  (it++ < maxit))
  {
    /* find zero crossing of line segment [(s0,G0), (s1,G1)] */
    double root = sqrt(G2*G2-G0*G1);
    if (fabs(root)<EPS)
      return -1.0; /* Hmmm */

    double sgn  = G0>G1 ? 1.0 : -1.0;
    s3          = s2 + ( s2-s0 )*sgn*G2/root;
    G3          = G(s3, data);


    /* if     G2*G3<0  */
    if (G2*G3 <= 0.0)
    {
      if (G2 > G3)
      {
        s0 = s3;
        G0 = G3;
        s1 = s2;
        G1 = G2;
      }
      else
      {
        s0 = s2;
        G0 = G2;
        s1 = s3;
        G1 = G3;

      }

    }
    else if (G0*G3 <= 0.0)
    {
      s1 = s3;
      G1 = G3;
    }
    else if (G1*G3 <= 0.0)
    {
      s0 = s3;
      G0 = G3;
    }
    else
    {
      print(
            "In ridder:\n"
            "G0=%10.10f, G1=%10.10f, G3=%10.10f\n",
            G0, G1, G3
            );
      getchar();
    }
    s2   = 0.5*(s0+s1);
    G2   = G(s2, data);
  }

  *iterations = it;
  return s3;
}



/* Start with bracket [0,1] with G(0)*G(1)<0.  Search by finding zero crossing
   sN of line segment [(sL,GL), (sR,GR)].  Set SL=sN if G(sN<0), sR=sN
   otherwise.*/
/*---------------------------------------------------------------------------*/
double
regulafalsi (double (*G)(double, void*), void *data, double tol, int maxit,
             int *iterations)
/*---------------------------------------------------------------------------*/
{
  *iterations = 0;

  double sn = ((double *) data)[2]; /* "Undefined" value */
  double Gn = G(sn, data);
  if (fabs(Gn) < tol) return sn;

  /* Initial guess is interval [0,1] of course */
  double s0=0.0;
  double G0=G(s0, data);
  if (fabs(G0)<EPS) return s0;

  double s1=1.0;
  double G1=G(s1, data);
  if (fabs(G1)<EPS) return s1;

  if (G0*G1 > 0) {
    print(no_root_str, "regulafalsi", G0, G1);
    return -1.0;
  }

  if (G0>G1) {
    double swap;
    swap = G0;
    G0   = G1;
    G1   = swap;
  }

  int it = 0;
  while ( (fabs(Gn) > tol)  &&  (it++ < maxit))  {

#if 0
    /* Unmodified Regula-Falsi */
    /* maintain bracket with G1>G0 */
    if ( Gn>0 ) {
      G1 = Gn;
      s1 = sn;
    }
    else {
      G0 = Gn;
      s0 = sn;
    }
#else
    /* Modified Regula-Falsi*/
    if ((Gn>0.0)==(G0>0.0)) {
      s0 = s1;
      G0 = G1;
    }
    else {
      /* const double gamma_illinois = 0.5; */
      const double gamma_pegasus  = G1/(G1+Gn);
      G0 *= gamma_pegasus;
    }
    s1 = sn;
    G1 = Gn;
#endif

    sn = s0 - (s1-s0)*G0/(G1-G0);
    Gn = G(sn, data);
  }

  *iterations = it;
  return sn;
}



/* Start with bracket [0,1] with G(0)*G(1)<0.  Search by finding sN=0.5(sL+sR).
   Set SL=sN if G(sN<0), sR=sN otherwise.*/
/*---------------------------------------------------------------------------*/
double
bisection (double (*G)(double, void*), void *data, double tol, int maxit,
           int *iterations)
/*---------------------------------------------------------------------------*/
{
  *iterations = 0;

  double sn = ((double *) data)[2];
  double Gn = G(sn, data);
  if (fabs(Gn)<tol) return sn;

  /* Initial guess is interval [0,1] of course */
  double s0=0.0;
  double G0=G(s0, data);
  if (fabs(G0)<EPS) return s0;

  double s1=1.0;
  double G1=G(s1, data);
  if (fabs(G1)<EPS) return s1;

  if (G0*G1 > 0.0) {
    print(no_root_str, "bisection", G0, G1);
    return -1.0;
  }

  if (G0>G1) {
    double swap;
    swap = G0;
    G0   = G1;
    G1   = swap;
  }

  int it=0;
  while ( (fabs(Gn)>tol) && (it++ < maxit) ) {
    if ( Gn>0 ) {
      G1 = Gn;
      s1 = sn;
    }
    else {
      G0 = Gn;
      s0 = sn;
    }

    sn = 0.5*(s0+s1);
    Gn = G(sn, data);
  }
  *iterations = it;
  if (it >= maxit) print("Warning: convergence criterion not met\n");
  return sn;
}

/* Local Variables:    */
/* c-basic-offset:4    */
/* End:                */
