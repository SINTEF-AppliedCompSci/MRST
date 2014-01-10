/*
Copyright (C) 2012 (c) Jostein R. Natvig <jostein natvig at gmail.com>

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is furnished to do
so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#ifdef MATLAB_MEX_FILE
#include "nlsolvers.h"

#include <mex.h>
extern int interrupt_signal;
#define print   mexPrintf
#define malloc  mxMalloc
#define calloc  mxCalloc
#define realloc mxRealloc
#define free    mxFree
#else

#define print printf
#include "nlsolvers.h"

#endif

static const char no_root_str[]=
    "  In %s:\n"
    "  With G(0) =% 5f, G(1) =% 5f, G(s) cannot have a zero in [0,1]!\n";



/*---------------------------------------------------------------------------*/
double
find_zero (double (*f)(double, void*), void *data, struct NonlinearSolverCtrl *ctrl)
/*---------------------------------------------------------------------------*/
{
    double zero;
    switch (ctrl->method) {
    default:
    case BISECTION:
        zero = bisection  (f, data, ctrl);
        break;

    case RIDDERS:
        zero = ridders    (f, data, ctrl);
        break;

    case REGULAFALSI:
        zero = regulafalsi(f, data, ctrl);
        break;
    }

    return zero;
}

/* Start with bracket [0,1] with G(0)*G(1)<0.  */
/*---------------------------------------------------------------------------*/
double
ridders (double (*G)(double, void*), void *data, struct NonlinearSolverCtrl *ctrl)
/*---------------------------------------------------------------------------*/
{
    double G0, G1, G2, G3;
    double s0, s1, s2, s3;
    double swap, sgn, root;
    int it;

    ctrl->iterations = 0;

    s2 = ctrl->initialguess;
    G2 = G(s2, data);
    if (fabs(G2) < ctrl->nltolerance) { return s2; }

    s0 = 0.0;
    G0 = G(s0, data);
    if (fabs(G0) < ctrl->nltolerance) { return s0; }

    s1 = 1.0;
    G1 = G(s1, data);
    if (fabs(G1) < ctrl->nltolerance) { return s1; }

    if (G0*G1 > 0)
    {
        print(no_root_str, "ridder", G0, G1);
        return -1.0;
    }

    if (G0>G1)
    {
        swap = G0;
        G0   = G1;
        G1   = swap;
    }

    s3 = 0;
    G3 = 10;

    it = 0;
    while ( (fabs(G3) > ctrl->nltolerance)  &&
            (ctrl->iterations++ < ctrl->maxiterations))
    {
        /* find zero crossing of line segment [(s0,G0), (s1,G1)] */
        root = sqrt(G2*G2 - G0*G1);

        sgn = G0>G1 ? 1.0 : -1.0;
        s3  = s2 + ( s2-s0 )*sgn*G2/root;
        G3  = G(s3, data);


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
            print("In ridder:\nG0=%10.10f, G1=%10.10f, "
                  "G3=%10.10f\n", G0, G1, G3);
        }
        s2   = 0.5*(s0+s1);
        G2   = G(s2, data);
    }

    ctrl->residual = G3;
    return s3;
}



/* Start with bracket [0,1] with G(0)*G(1)<0.  Search by finding zero crossing
   sN of line segment [(sL,GL), (sR,GR)].  Set SL=sN if G(sN<0), sR=sN
   otherwise.*/
/*---------------------------------------------------------------------------*/
double
regulafalsi (double (*G)(double, void*), void *data, struct NonlinearSolverCtrl *ctrl)
/*---------------------------------------------------------------------------*/
{
    double Gn, G0, G1;
    double sn, s0, s1;
    double swap, gamma_pegasus;
    int it;

    ctrl->iterations = 0;

    sn = ctrl->initialguess;
    Gn = G(sn, data);
    if (fabs(Gn) < ctrl->nltolerance) { return sn; }

    /* Initial guess is interval [0,1] of course */
    s0 = 0.0;
    G0 = G(s0, data);
    if (fabs(G0) < ctrl->nltolerance) { return s0; }

    s1 = 1.0;
    G1 = G(s1, data);
    if (fabs(G1) < ctrl->nltolerance) { return s1; }

    if (G0*G1 > 0)
    {
        print(no_root_str, "regulafalsi", G0, G1);
        return -1.0;
    }

    if (G0>G1)
    {
        swap = G0;
        G0   = G1;
        G1   = swap;
    }

    it = 0;
    while ( (fabs(Gn) > ctrl->nltolerance)  &&
            (ctrl->iterations++ < ctrl->maxiterations))
    {

#if 0
        /* Unmodified Regula-Falsi */
        /* maintain bracket with G1>G0 */
        if ( Gn>0 )
        {
            G1 = Gn;
            s1 = sn;
        }
        else
        {
            G0 = Gn;
            s0 = sn;
        }
#else
        /* Modified Regula-Falsi*/
        if ((Gn>0.0)==(G0>0.0))
        {
            s0 = s1;
            G0 = G1;
        }
        else
        {
            /* const double gamma_illinois = 0.5; */
            gamma_pegasus  = G1/(G1+Gn);
            G0 *= gamma_pegasus;
        }
        s1 = sn;
        G1 = Gn;
#endif

        sn = s0 - (s1-s0)*G0/(G1-G0);
        Gn = G(sn, data);
    }

    ctrl->residual = Gn;
    return sn;
}



/* Start with bracket [0,1] with G(0)*G(1)<0.  Search by finding sN=0.5(sL+sR).
   Set SL=sN if G(sN<0), sR=sN otherwise.*/
/*---------------------------------------------------------------------------*/
double
bisection (double (*G)(double, void*), void *data, struct NonlinearSolverCtrl *ctrl)
/*---------------------------------------------------------------------------*/
{
    double Gn, G0, G1;
    double sn, s0, s1;
    double swap;
    int it;

    ctrl->iterations = 0;

    sn = ctrl->initialguess;
    Gn = G(sn, data);
    if (fabs(Gn) < ctrl->nltolerance) { return sn; }

    /* Initial guess is interval [0,1] of course */
    s0 = 0.0;
    G0 = G(s0, data);
    if (fabs(G0) < ctrl->nltolerance) { return s0; }

    s1 = 1.0;
    G1 = G(s1, data);
    if (fabs(G1) < ctrl->nltolerance) { return s1; }

    if (G0*G1 > 0.0)
    {
        print(no_root_str, "bisection", G0, G1);
        return -1.0;
    }

    if (G0>G1)
    {
        swap = G0;
        G0   = G1;
        G1   = swap;
    }

    it=0;
    while ( (fabs(Gn)>ctrl->nltolerance) &&
            (ctrl->iterations++ < ctrl->maxiterations) )
    {
        if ( Gn>0 )
        {
            G1 = Gn;
            s1 = sn;
        }
        else
        {
            G0 = Gn;
            s0 = sn;
        }

        sn = 0.5*(s0+s1);
        Gn = G(sn, data);
    }
    if (ctrl->iterations >= ctrl->maxiterations)
    {
        print("Warning: convergence criterion not met\n");
    }
    ctrl->residual = Gn;
    return sn;
}


/* Local Variables:    */
/* c-basic-offset:4    */
/* End:                */
