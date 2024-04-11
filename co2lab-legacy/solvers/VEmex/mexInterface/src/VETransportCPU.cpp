/*
-------------------------------------------------------------------------------
File Name : explicitTransport.cpp

Purpose :

Creation Date : 2010-10-29

Last Modified : fr. 03. des. 2010 kl. 17.56 +0100

Created By :  Martin Ertsaas (martiert@student.matnat.uio.no)
-------------------------------------------------------------------------------
*/

// Hack to compile on windows with Visual Studio
#if (_MSC_VER >= 1600)
#include <yvals.h>
#define __STDC_UTF_16__
#endif

#include "mex.h"

#include <string.h>

#include "extractStructs.h"
#include "State.h"
#include "Grid2D.h"
#include "Rock.h"
#include "Fluid.h"
#include "Options.h"
#include "VESimulator.h"

// Another VS hack
#if (_MSC_VER >= 1600)
#pragma comment(lib, "libmx.lib")
#pragma comment(lib, "libmat.lib")
#pragma comment(lib, "libmex.lib")
#pragma comment(lib, "libmatlb.lib")
#endif

static VESimulator<double> *simud = NULL;
static int initd = 0;
static VESimulator<float> *simuf = NULL;
static int initf = 0;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if(nrhs == 0) {
        if(simud) {
            delete simud;
            simud = NULL;
        }
        initd = 0;

        if(simuf) {
            delete simuf;
            simuf = NULL;
        }
        initf = 0;

        return;
    }

    if(nrhs < 5)
        mexErrMsgTxt("Function takes at least 5 input arguments");

    if(!mxIsStruct(prhs[0]))
        mexErrMsgTxt("First input arguments needs to be a 2D state struct.");

    if(!mxIsStruct(prhs[1]))
        mexErrMsgTxt("Second input argument needs to be a 2D grid struct.");

    if(!mxIsStruct(prhs[3]))
        mexErrMsgTxt("Fourth input argument needs to be a 2D rock struct.");

    if(!mxIsStruct(prhs[4]))
        mexErrMsgTxt("Fifth argument needs to be a 2D fluid struct.");

    if((nrhs - 5)%2 != 0)
        mexErrMsgTxt("The optional arguments needs to come in 'key'/value \
                pairs.\n");


    opt_t *options  = extractOptions(prhs, nrhs);

    omp_set_num_threads(options->flag);
    //omp_set_dynamic(0);//seamed to make things run slower

    if(options->flag != 0) {
        State<double> *state    = extractStateStruct<double>(prhs[0]);
        double dT = mxGetScalar(prhs[2]);

        if(options->dt == 0.0) {
            options->dt = dT;
        }

        if(!initd) {
            Grid2D<double> *grid = extractGrid2DStruct<double>(prhs[1]);
            Rock<double> *rock   = extractRockStruct<double>(prhs[3]);
            Fluid<double> *fluid = extractFluidStruct<double>(prhs[4]);
            simud = new VESimulator<double>(state, grid, rock, fluid, options);
            initd = 1;
        } else {
            simud->changeState(state);
            simud->changeOptions(options);
        }

        simud->doTimeStep(dT);

        state = simud->getState();

        mxArray *height = mxCreateDoubleMatrix(state->numCells, 1, mxREAL);
        mxArray *maxheight = mxCreateDoubleMatrix(state->numCells, 1, mxREAL);
        double *h = mxGetPr(height);
        double *h_max = mxGetPr(maxheight);

        state->getHeight(h, h_max);

        plhs[0] = height;
        plhs[1] = maxheight;
    } else {
        omp_set_num_threads(omp_get_max_threads());

        State<float> *state    = extractStateStruct<float>(prhs[0]);

        float dT       = mxGetScalar(prhs[2]);

        if(options->dt == 0.0) {
            options->dt = dT;
        }

        if(!initf) {
            Grid2D<float> *grid = extractGrid2DStruct<float>(prhs[1]);
            Rock<float> *rock = extractRockStruct<float>(prhs[3]);
            Fluid<float> *fluid = extractFluidStruct<float>(prhs[4]);
            simuf = new VESimulator<float>(state, grid, rock, fluid, options);
            initf = 1;
        } else {
            simuf->changeState(state);
            simuf->changeOptions(options);
        }

        simuf->doTimeStep(dT);

        state = simuf->getState();

        mxArray *height = mxCreateDoubleMatrix(state->numCells, 1, mxREAL);
        mxArray *maxheight = mxCreateDoubleMatrix(state->numCells, 1, mxREAL);
        double *h = mxGetPr(height);
        double *h_max = mxGetPr(maxheight);

        state->getHeight(h, h_max);

        plhs[0] = height;
        plhs[1] = maxheight;
    }
    free(options);
}
