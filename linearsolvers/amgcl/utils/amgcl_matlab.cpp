//
// include necessary system headers
//
#include <mex.h>
#include "matrix.h"
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
    const mxArray * pa;
    
    if (nrhs != 3) { 
	    mexErrMsgTxt("3 input arguments required."); 
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

    pa = prhs[2];
    double tolerance = mxGetScalar(mxGetField(pa, 0, "tolerance"));
    int maxiter = mxGetScalar(mxGetField(pa, 0, "maxIterations"));
    int coarsen_id = mxGetScalar(mxGetField(pa, 0, "coarsening"));
    int relax_id = mxGetScalar(mxGetField(pa, 0, "relaxation"));
    int solver_id = mxGetScalar(mxGetField(pa, 0, "solver"));
    int precond_id = mxGetScalar(mxGetField(pa, 0, "preconditioner"));
    
    std::vector<double> b(n);
    for(int ix = 0; ix < n; ix++){
        b[ix] = rhs[ix];
    }
    int M = (int)m;
    /***************************************
     * Start AMGCL-link and select options *
     ***************************************/
    typedef amgcl::backend::builtin<double> Backend;
    
    boost::property_tree::ptree prm;
    /* Select preconditioning strategy */
    switch(precond_id) {
        case 1:
            relaxParam = "precond.relax.type";
            prm.put("precond.class", amgcl::runtime::precond_class::amg);
            break;
        case 2: 
            relaxParam = "precond.type";
            prm.put("precond.class", amgcl::runtime::precond_class::relaxation);
            break;
        case 3:
            relaxParam = "precond.relax.type";
            prm.put("precond.class", amgcl::runtime::precond_class::dummy);
            break;
        default : mexErrMsgTxt("Unknown precond_id."); 
    }
    /* AMG specific options */
    if(precond_id == 1){
        /* Select coarsening strategy */
        switch(coarsen_id) {
            case 1: 
                prm.put("precond.coarsening.type",  amgcl::runtime::coarsening::smoothed_aggregation);
                break;
            case 2: 
                prm.put("precond.coarsening.type",  amgcl::runtime::coarsening::ruge_stuben);
                break;
            case 3: 
                prm.put("precond.coarsening.type",  amgcl::runtime::coarsening::aggregation);
                break;
            case 4: 
                prm.put("precond.coarsening.type",  amgcl::runtime::coarsening::smoothed_aggr_emin);
                break;
            default : mexErrMsgTxt("Unknown coarsen_id."); 
        }
        /* When is a level coarse enough */
        int coarse_enough = mxGetScalar(mxGetField(pa, 0, "coarse_enough"));
        if (coarse_enough >= 0){
            prm.put("precond.coarse_enough", coarse_enough);
        }
        /* Use direct solver for coarse sys */
        bool direct_coarse = mxGetScalar(mxGetField(pa, 0, "direct_coarse"));
        prm.put("precond.direct_coarse", direct_coarse);
        /* Max levels */
        int max_levels = mxGetScalar(mxGetField(pa, 0, "max_levels"));
        if (max_levels >= 0){
            prm.put("precond.max_levels", max_levels);
        }
        /* Number of cycles */
        int ncycle = mxGetScalar(mxGetField(pa, 0, "ncycle"));
        if (ncycle >= 0){
            prm.put("precond.ncycle", ncycle);
        }
        /* Pre cycles */
        int npre = mxGetScalar(mxGetField(pa, 0, "npre"));
        if (npre >= 0){
            prm.put("precond.npre", npre);
        }
        /* Post cycles */
        int npost = mxGetScalar(mxGetField(pa, 0, "npost"));
        if (npost >= 0){
            prm.put("precond.npost", npost);
        }
        /* Pre cycles (precond) */
        int pre_cycles = mxGetScalar(mxGetField(pa, 0, "pre_cycles"));
        if (pre_cycles >= 0){
            prm.put("precond.pre_cycles", pre_cycles);
        }
    }
    /* Select relaxation strategy */
    switch(relax_id) {
        case 1: 
            prm.put(relaxParam,  amgcl::runtime::relaxation::spai0);
            break;
        case 2: 
            prm.put(relaxParam,  amgcl::runtime::relaxation::gauss_seidel);
            break;
        case 3: 
            prm.put(relaxParam,  amgcl::runtime::relaxation::ilu0);
            break;
        case 4: 
            prm.put(relaxParam,  amgcl::runtime::relaxation::iluk);
            break;
        case 5: 
            prm.put(relaxParam,  amgcl::runtime::relaxation::ilut);
            break;
        case 6: 
            prm.put(relaxParam,  amgcl::runtime::relaxation::damped_jacobi);
            break;
        case 7: 
            prm.put(relaxParam,  amgcl::runtime::relaxation::spai1);
            break;
        case 8: 
            prm.put(relaxParam,  amgcl::runtime::relaxation::chebyshev);
            break;
        default : mexErrMsgTxt("Unknown relax_id."); 
    }

    /* Select solver */
    switch(solver_id) {
        case 1: 
            prm.put("solver.type",  amgcl::runtime::solver::bicgstab);
            break;
        case 2: 
            prm.put("solver.type",  amgcl::runtime::solver::cg);
            break;
        case 3: 
            prm.put("solver.type",  amgcl::runtime::solver::bicgstabl);
            break;
        case 4: 
            prm.put("solver.type",  amgcl::runtime::solver::gmres);
            break;
        case 5: 
            prm.put("solver.type",  amgcl::runtime::solver::lgmres);
            break;
        case 6: 
            prm.put("solver.type",  amgcl::runtime::solver::fgmres);
            break;
        case 7: 
            prm.put("solver.type",  amgcl::runtime::solver::idrs);
            break;
        default : mexErrMsgTxt("Unknown solver_id."); 
    }
    /* Set tolerance */
    prm.put("solver.tol", tolerance);
    if(maxiter > 0){
        prm.put("solver.maxiter", maxiter);
    }
    // TODO: Expose more parameters through generic cell-array style interface?
    
    amgcl::make_solver<
        amgcl::runtime::preconditioner<Backend>,
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

