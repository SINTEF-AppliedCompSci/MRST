//
// include necessary system headers
//
#define _USE_MATH_DEFINES
#include "matrix.h"
#include <array>
#include <cmath>
#include <mex.h>


#include <fstream>
#include <iostream>
/* MEX interfaces */
/* Block system support */
#include "FlexibleSolver.hpp"
#include <cmath>
#include <iomanip>
#include <limits>
//#include <boost/program_options.hpp>
// #include <dune/common/parallel/mpihelper.hh>
// #include <dune/istl/bcrsmatrix.hh>
// #include <dune/istl/matrixmarket.hh>
// #include <dune/common/fmatrix.hh>
// #include <dune/istl/solvers.hh>
// #include <dune/istl/preconditioners.hh>
// #include <dune/istl/umfpack.hh>
// #include <dune/istl/solvers.hh>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>

void
mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])

{
    char* matrixfilename = mxArrayToString(prhs[0]);
    char* rhsfilename = mxArrayToString(prhs[1]);
    // int m  = mxArrayToInt(prhs[2]);
    double mm = mxGetScalar(prhs[2]);
    int m(mm);
    double bzm = mxGetScalar(prhs[3]);
    int bz(bzm);
    /* Assign a pointer to the output */
    // y = mxGetDoubles(plhs[0]);
    std::string rhsfile(rhsfilename);
    std::string matrixfile(matrixfilename);
    plhs[0] = mxCreateDoubleMatrix(m, 1, mxREAL);
    double* result = mxGetPr(plhs[0]);
    namespace pt = boost::property_tree;
    pt::ptree prm;
    {
        char* chopt = mxArrayToString(prhs[4]);
        std::istringstream str(chopt);
        pt::read_json(str, prm);
        std::cout << "options in prm" << std::endl;
        std::ofstream file("options.json");
        // pt::write_json(std::cout, prm);
        pt::write_json(file, prm);
    }
    if (bz == 1) {
        Dune::FlexibleSolver<1> solver(prm);
        solver.solve(result, matrixfile, rhsfile);
    } else if (bz == 2) {
        Dune::FlexibleSolver<2> solver(prm);
        solver.solve(result, matrixfile, rhsfile);
    } else if (bz == 3) {
        Dune::FlexibleSolver<3> solver(prm);
        solver.solve(result, matrixfile, rhsfile);
    } else {
        std::cout << "Flexible solver not implemented for blocksize " << bz << std::endl;
    }
    return;
}
