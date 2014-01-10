/*
  Copyright 2012 SINTEF ICT, Applied Mathematics.

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <mex.h>
#include "../mrst_api/grid.h"
#include "../mrst_api/mrst_api.h"
#include <opm/core/transport/reorder/TransportModelTwophase.hpp>
#include <opm/core/fluid/IncompPropertiesBasic.hpp>
#include <cmath>


void
mexFunction(int nlhs,       mxArray *plhs[],
            int nrhs, const mxArray *prhs[])
{
    if (nrhs == 6 && nlhs == 1) {
        UnstructuredGrid* g = mrst_grid  (prhs[0]);
        double* porevol     = mxGetPr    (prhs[1]);
        double* flux        = mxGetPr    (prhs[2]);
        double* src         = mxGetPr    (prhs[3]);
        double dt           = mxGetScalar(prhs[4]);
        double* s0          = mxGetPr    (prhs[5]);

        Opm::parameter::ParameterGroup param;
        param.disableOutput();
        Opm::IncompPropertiesBasic props(param, g->dimensions, g->number_of_cells);
        const double tol = 1e-8;
        const int maxit = 30;
        Opm::TransportModelTwophase transport(*g, props, tol, maxit);

        // Since the solver overwrites saturation, we must take care to
        // preserve it properly. Also, we send
        std::vector<double> sat(s0, s0 + 2*g->number_of_cells);
        if (std::fabs(sat[0] + sat[1] - 1.0) > 1e-10) {
            mexErrMsgTxt("Saturations in first cell do not sum to 1.\n");
            return;
        }
        transport.solve(flux, porevol, src, dt, sat);
        plhs[0] = mxCreateDoubleMatrix(2, g->number_of_cells, mxREAL);
        double* sat_out = mxGetPr(plhs[0]);
        std::copy(sat.begin(), sat.end(), sat_out);
        free_mrst_grid(g);
    } else {
        mexErrMsgTxt("Wrong number of inputs (need 6) or outputs (need 1).\n");
    }
}
