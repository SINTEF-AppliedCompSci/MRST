//
// include necessary system headers
//
#include <mex.h>
#include <iostream>
#include <amgcl/make_solver.hpp>
#include <amgcl/solver/bicgstab.hpp>
#include <amgcl/amg.hpp>
#include <amgcl/coarsening/smoothed_aggregation.hpp>
#include <amgcl/relaxation/spai0.hpp>
#include <amgcl/relaxation/runtime.hpp>


#include <amgcl/adapter/crs_tuple.hpp>
#include <amgcl/adapter/zero_copy.hpp>

#include <amgcl/backend/builtin.hpp>
#include <amgcl/runtime.hpp>

#include <amgcl/preconditioner/dummy.hpp>
#include <amgcl/preconditioner/runtime.hpp>




#include <string>

#include <boost/program_options.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/foreach.hpp>

#include <amgcl/preconditioner/cpr.hpp>
#include <amgcl/adapter/crs_tuple.hpp>
//#include <amgcl/io/mm.hpp>
//#include <amgcl/io/binary.hpp>
//#include <amgcl/profiler.hpp>

void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray *prhs[] )
     
{ 
    double *result; 
    double *rhs;
    double *err;
    mwSize m,n,nnz;
    mwIndex * cols;
    mwIndex * rows;
    double * entries;
    std::string relaxParam;
    
    if (nrhs != 4) { 
	    mexErrMsgTxt("4 input arguments required."); 
    } else if (nlhs > 2) {
	    mexErrMsgTxt("Wrong number of output arguments."); 
    } 

    m = mxGetM(prhs[0]); 
    n = mxGetN(prhs[0]);
    if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) ||  !mxIsSparse(prhs[0]) ) { 
	    mexErrMsgTxt("Matrix should be a real sparse matrix."); 
        return;
    } 
    if (!mxIsDouble(prhs[1]) || mxIsComplex(prhs[1])) { 
	    mexErrMsgTxt("Right hand side must be real double column vector.");
        return; 
    } 
    // main();
    plhs[0] = mxCreateDoubleMatrix(m, 1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
    
    result = mxGetPr(plhs[0]);
    err = mxGetPr(plhs[1]);
    
    cols    = mxGetJc(prhs[0]);
    rows    = mxGetIr(prhs[0]);
    entries = mxGetPr(prhs[0]);
    nnz  = mxGetNzmax(prhs[0]);
    rhs     = mxGetPr(prhs[1]);

    double tolerance = mxGetScalar(prhs[2]);
    int maxiter = mxGetScalar(prhs[3]);

    std::vector<double> b(n);
    for(int ix = 0; ix < n; ix++){
        b[ix] = rhs[ix];
    }
    int M = (int)m;
    /***************************************
     * Start AMGCL-link and select options *
     ***************************************/
    typedef amgcl::backend::builtin<double> Backend;
    
    typedef
        amgcl::runtime::amg<Backend>
        PPrecond;

    typedef
        amgcl::runtime::relaxation::as_preconditioner<Backend>
        SPrecond;


    
    boost::property_tree::ptree prm;
    /* Set tolerance */
    prm.put("solver.tol", tolerance);
    if(maxiter > 0){
        prm.put("solver.maxiter", maxiter);
    }
    prm.put("precond.block_size", 2);

    amgcl::make_solver<
        amgcl::preconditioner::cpr<PPrecond, SPrecond>,
        amgcl::runtime::iterative_solver<Backend>
        > solve(amgcl::adapter::zero_copy(n, &cols[0], &rows[0], &entries[0]), prm);
    /***************************************
     *        Solve problem                *
     ***************************************/
    std::vector<double> x(M, 0.0);
    int    iters;
    double error;
    boost::tie(iters, error) = solve(b, x);
    
    for(int ix=0; ix < M; ix++){
        result[ix] = x[ix];
    }
    err[0] = error;

    return;
}

