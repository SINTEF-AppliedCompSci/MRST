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
#include <opm/core/transport/reorder/TransportModelTracerTof.hpp>
#include <opm/core/transport/reorder/TransportModelTracerTofDiscGal.hpp>
#include <opm/core/utility/ErrorMacros.hpp>
#include <string>
#include <vector>

#define MEX_MESSAGE(x) mexPrintf("\nError in file %s, line %s.\n", __FILE__, __LINE__)
#define MEX_THROW(x) do { MEX_MESSAGE(x); throw std::exception(); } while(0)

// Note: this function has not been tested properly with disc_method == DG0 or DG1,
// since it triggers LAPACK calls, and I have been unable to prevent crashes due
// to Matlab using its own version of BLAS and LAPACK.
// Calling this function with 4 input arguments has been tested.

void
mexFunction(int nlhs,       mxArray *plhs[],
            int nrhs, const mxArray *prhs[])
{
    if ((nrhs == 4 || nrhs == 5 || nrhs == 6) && nlhs == 1) {
        UnstructuredGrid* g = mrst_grid  (prhs[0]);
        double* porevol     = mxGetPr    (prhs[1]);
        double* flux        = mxGetPr    (prhs[2]);
        double* src         = mxGetPr    (prhs[3]);

        // Get discretization method.
        enum DiscMethod { FV, DG0, DG1 };
        DiscMethod disc_method = FV;
        if (nrhs >= 5) {
            const int MaxStrLen = 100;
            char buf[MaxStrLen];
            int notok = mxGetString(prhs[4], buf, MaxStrLen);
            if (notok) {
                MEX_THROW("Error in string conversion.");
            }
            std::string dm(buf);
            if (dm == "FV") {
                disc_method = FV;
            } else if (dm == "DG0") {
                disc_method = DG0;
            } else if (dm == "DG1") {
                disc_method = DG1;
            } else {
                MEX_THROW("Unknown discretization method " << dm);
            }
        }

        // Get velocity interpolation method.
        enum VelInterpMethod { Constant, ECVI };
        VelInterpMethod vi_method = Constant;
        if (nrhs >= 6) {
            const int MaxStrLen = 100;
            char buf[MaxStrLen];
            int notok = mxGetString(prhs[5], buf, MaxStrLen);
            if (notok) {
                MEX_THROW("Error in string conversion.");
            }
            std::string vm(buf);
            if (vm == "Constant") {
                vi_method = Constant;
            } else if (vm == "ECVI") {
                vi_method = ECVI;
            } else {
                MEX_THROW("Unknown velocity interpolation method " << vm);
            }
        }

        // Compute TOF.
        std::vector<double> tof;
        int bf_per_cell = -1;
        switch (disc_method) {
        case FV:
            {
                bf_per_cell = 1;
                Opm::TransportModelTracerTof tof_computer(*g);
                tof_computer.solveTof(flux, porevol, src, tof);
                break;
            }
        case DG0:
            {
                bf_per_cell = 1;
                Opm::TransportModelTracerTofDiscGal tof_computer(*g, false);
                tof_computer.solveTof(flux, porevol, src, 0, tof);
                break;
            }
        case DG1:
            {
                bf_per_cell = 4;
                Opm::TransportModelTracerTofDiscGal tof_computer(*g, vi_method == ECVI);
                tof_computer.solveTof(flux, porevol, src, 1, tof);
                break;
            }
        default:
            MEX_THROW("Unknown discretization method index " << disc_method);
        }

        // Write output.
        if (tof.size() != bf_per_cell*g->number_of_cells) {
            MEX_THROW("Error: tof data size not equal to " << bf_per_cell << " per cell.");
        }
        plhs[0] = mxCreateDoubleMatrix(bf_per_cell, g->number_of_cells, mxREAL);
        double* tof_out = mxGetPr(plhs[0]);
        std::copy(tof.begin(), tof.end(), tof_out);
        free_mrst_grid(g);
    } else {
        mexErrMsgTxt("Wrong number of inputs (need 4) or outputs (need 1).\n");
    }
}
