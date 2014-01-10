#include <math.h>

#ifdef MATLAB_MEX_FILE
#include <mex.h>
#include "fluid.h"
#else
#include <stdio.h>
#include <fluid.h>
#endif


static int    nr;
static double *viscw;
static double *visco;
static double *nw;
static double *no;
static double *srw;
static double *sro;

/*---------------------------------------------------------------------------*/
double fluxfun(double sw, const int reg)
/*---------------------------------------------------------------------------*/
{
    double den, so, mw, mo;

#ifdef MATLAB_MEX_FILE
    mxAssert ((0 <= reg) && (reg < nr),"Region number out of bounds.");
#else
    if(0 <= reg)
    {
        fprintf(stderr, "Region number out of bounds.");
        exit(EXIT_FAILURE);
    }
#endif
    den = 1 - srw[reg] - sro[reg];
    so  = (1 - sw - sro[reg])/den;
    sw  = (sw - srw[reg])/den;
    sw  = sw > 1.0 ? 1.0 : sw;
    sw  = sw < 0.0 ? 0.0 : sw;

    mw  = pow(sw, nw[reg]) / viscw[0];
    mo  = pow(so, no[reg]) / visco[0];

    return mw / (mw + mo);
}

/*---------------------------------------------------------------------------*/
double dfluxfun(double sw, const int reg)
/*---------------------------------------------------------------------------*/
{
    double den, so, mw, mo, dmw, dmo, fw;

    den = 1 - srw[reg] - sro[reg];
    so = (1 - sw - sro[reg])/den;
    sw = (sw - srw[reg])/den;
    sw = sw > 1.0 ? 1.0 : sw;
    sw = sw < 0.0 ? 0.0 : sw;

    mw  = pow(sw, nw[reg]) / viscw[0];
    mo  = pow(so, no[reg]) / visco[0];

    dmw = nw[reg] * pow(sw, nw[reg]-1) / viscw[0];
    dmo = no[reg] * pow(so, no[reg]-1) / visco[0];

    fw = mw / (mw+mo);

    return (dmw * (1-fw) - dmo * fw) / (mw+mo);
}

#ifdef MATLAB_MEX_FILE
/* ------------------------------------------------------------------ */
static const mxArray*
getField(const mxArray *a, const char *field)
/* ------------------------------------------------------------------ */
{
    int fld_no = mxGetFieldNumber(a, field);

    if (fld_no < 0)
    {
       mexErrMsgTxt("Missing field in fluid opts.");
    }
    return mxGetFieldByNumber(a , 0, fld_no);
}


/*---------------------------------------------------------------------------*/
void init_fluid(const mxArray *arr)
/*---------------------------------------------------------------------------*/
{
    viscw = mxGetPr(getField(arr, "viscw"));
    visco = mxGetPr(getField(arr, "visco"));
    srw   = mxGetPr(getField(arr, "srw"));
    sro   = mxGetPr(getField(arr, "sro"));
    nw    = mxGetPr(getField(arr, "nw"));
    no    = mxGetPr(getField(arr, "no"));

    nr    = mxGetNumberOfElements(getField(arr, "no"));
}
#else
void init_fluid(void)
{
    fprintf(stderr, "Not implemented");
}
#endif

/*---------------------------------------------------------------------------*/
int getNumSatRegions(void)
/*---------------------------------------------------------------------------*/
{
    return nr;
}

/* Local Variables:    */
/* c-basic-offset:4    */
/* End:                */
