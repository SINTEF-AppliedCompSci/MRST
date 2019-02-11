//
// include necessary system headers
//
#define _USE_MATH_DEFINES
#include <cmath>
#include <mex.h>
#include <array>
#include "matrix.h"


#include <iostream>
#include <amgcl/make_solver.hpp>
#include <amgcl/solver/runtime.hpp>

#include <amgcl/amg.hpp>
#include <amgcl/adapter/crs_tuple.hpp>
#include <amgcl/backend/block_crs.hpp>
#include <amgcl/adapter/zero_copy.hpp>
#include <amgcl/backend/builtin.hpp>
#include <amgcl/coarsening/runtime.hpp>
#include <amgcl/relaxation/runtime.hpp>
#include <amgcl/relaxation/as_preconditioner.hpp>
#include <string>
#include <chrono>
#include <amgcl/preconditioner/runtime.hpp>

#include <amgcl/preconditioner/cpr.hpp>
#include <amgcl/preconditioner/cpr_drs.hpp>

/* MEX interfaces */
#include "amgcl_mex_utils.cpp"

/* Block system support */
#ifndef AMGCL_BLOCK_SIZES
#  define AMGCL_BLOCK_SIZES (2)(3)(4)(5)
#endif

#ifndef SOLVER_BACKEND_BUILTIN
#  define SOLVER_BACKEND_BUILTIN
#endif
#include <amgcl/backend/builtin.hpp>
#include <amgcl/value_type/static_matrix.hpp>
#include <amgcl/adapter/block_matrix.hpp>
#include <amgcl/make_block_solver.hpp>
#include <amgcl/adapter/crs_tuple.hpp>

typedef amgcl::backend::builtin<double> Backend;

void solve_cpr(int n, mwIndex * cols, mwIndex * rows, double * entries, const mxArray * pa,
        std::vector<double> b, std::vector<double> & x, double tolerance,
        int maxiter, int & iters, double & error){

    int block_size = mxGetScalar(mxGetField(pa, 0, "block_size"));
    int active_rows = mxGetScalar(mxGetField(pa, 0, "active_rows"));
    bool use_drs = mxGetScalar(mxGetField(pa, 0, "use_drs"));


    int relax_p_id = mxGetScalar(mxGetField(pa, 0, "relaxation"));
    int relax_s_id = mxGetScalar(mxGetField(pa, 0, "s_relaxation"));

    bool verbose = mxGetScalar(mxGetField(pa, 0, "verbose"));


    /***************************************
     * Start AMGCL-link and select options *
     ***************************************/

    typedef
    amgcl::amg<Backend, amgcl::runtime::coarsening::wrapper, amgcl::runtime::relaxation::wrapper>
        PPrecond;

    typedef
    amgcl::relaxation::as_preconditioner<Backend, amgcl::runtime::relaxation::wrapper>
            SPrecond;
    boost::property_tree::ptree prm;
    /* Set tolerance */
    prm.put("solver.tol", tolerance);
    if(maxiter > 0){
        prm.put("solver.maxiter", maxiter);
    }
    prm.put("precond.block_size", block_size);
    prm.put("precond.active_rows", active_rows);
    /* Select coarsening strategy */
    amg_opts c_opt;
    setCoarseningStructMex(c_opt, pa);
    setCoarseningAMGCL(prm, "precond.pprecond.", c_opt);

    /* Select relaxation strategy for pressure solver */
    relax_opts pr_opt;
    // pr_opt.relax_id = relax_p_id;
    setRelaxationStructMex(pr_opt, pa, "");
    setRelaxationAMGCL(prm, "precond.pprecond.relax.", pr_opt);
    /* Select relaxation strategy for second stage solver */
    relax_opts ps_opt;
    // ps_opt.relax_id = relax_s_id;
    setRelaxationStructMex(ps_opt, pa, "s_");
    setRelaxationAMGCL(prm, "precond.sprecond.", ps_opt);

    /* Select solver */
    solver_opts sol_opt;
    setSolverStructMex(sol_opt, pa);
    setSolverAMGCL(prm, "solver.", sol_opt);

    /***************************************
     *        Solve problem                *
     ***************************************/
    auto t1 = std::chrono::high_resolution_clock::now();

    if(use_drs){
        double dd = mxGetScalar(mxGetField(pa, 0, "drs_eps_dd"));
        double ps = mxGetScalar(mxGetField(pa, 0, "drs_eps_ps"));
        prm.put("precond.eps_dd", dd);
        prm.put("precond.eps_ps", ps);

        mxArray * drs_weights_mx = mxGetField(pa, 0, "drs_row_weights");
        size_t drs_weights_n = mxGetM(drs_weights_mx);

        if(drs_weights_n>0){
            double * drs_weights = mxGetPr(drs_weights_mx);
            prm.put("precond.weights", drs_weights);
            prm.put("precond.weights_size", drs_weights_n);
        }
        amgcl::make_solver<
            amgcl::preconditioner::cpr_drs<PPrecond, SPrecond>,
            amgcl::runtime::solver::wrapper<Backend>
            > solve(amgcl::adapter::zero_copy(n, &cols[0], &rows[0], &entries[0]), prm);
        auto t2 = std::chrono::high_resolution_clock::now();

        if(verbose){
            std::cout << "CPR setup took "
                      << (double)std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count()/1000.0
                      << " seconds\n";
        }
        std::tie(iters, error) = solve(b, x);

        if(verbose){
            std::cout << solve << std::endl;
        }
    }else{
        amgcl::make_solver<
            amgcl::preconditioner::cpr<PPrecond, SPrecond>,
            amgcl::runtime::solver::wrapper<Backend>
            > solve(amgcl::adapter::zero_copy(n, &cols[0], &rows[0], &entries[0]), prm);

        auto t2 = std::chrono::high_resolution_clock::now();
        if(verbose){
            std::cout << "CPR setup took "
                      << (double)std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count()/1000.0
                      << " seconds\n";
        }
        std::tie(iters, error) = solve(b, x);

        if(verbose){
            std::cout << solve << std::endl;
        }
    }
}

void solve_regular(int n, mwIndex * cols, mwIndex * rows, double * entries, const mxArray * pa,
        std::vector<double> b, std::vector<double> & x, double tolerance,
        int maxiter, int & iters, double & error){

    int relax_id = mxGetScalar(mxGetField(pa, 0, "relaxation"));
    bool verbose = mxGetScalar(mxGetField(pa, 0, "verbose"));
    int precond_id = mxGetScalar(mxGetField(pa, 0, "preconditioner"));
    std::string relaxParam;
    /***************************************
     * Start AMGCL-link and select options *
     ***************************************/
    typedef amgcl::backend::builtin<double> Backend;

    typedef
    amgcl::amg<Backend, amgcl::runtime::coarsening::wrapper, amgcl::runtime::relaxation::wrapper>
        PPrecond;

    typedef
    amgcl::relaxation::as_preconditioner<Backend, amgcl::runtime::relaxation::wrapper>
            SPrecond;

    boost::property_tree::ptree prm;
    /* Set tolerance */
    prm.put("solver.tol", tolerance);
    if(maxiter > 0){
        prm.put("solver.maxiter", maxiter);
    }
    switch(precond_id) {
        case 1:
            relaxParam = "precond.relax.";
            prm.put("precond.class", amgcl::runtime::precond_class::amg);
            break;
        case 2:
            relaxParam = "precond.";
            prm.put("precond.class", amgcl::runtime::precond_class::relaxation);
            break;
        case 3:
            relaxParam = "precond.relax.";
            prm.put("precond.class", amgcl::runtime::precond_class::dummy);
            break;
        default : mexErrMsgTxt("Unknown precond_id.");
    }

    if(precond_id == 1){
        /* Select coarsening strategy */
        amg_opts c_opt;
        setCoarseningStructMex(c_opt, pa);
        setCoarseningAMGCL(prm, "precond.", c_opt);
    }
    /* Select relaxation strategy for solver */
    relax_opts pr_opt;
    // pr_opt.relax_id = relax_id;
    setRelaxationStructMex(pr_opt, pa, "");
    setRelaxationAMGCL(prm, relaxParam, pr_opt);

    /* Select solver */
    solver_opts sol_opt;
    setSolverStructMex(sol_opt, pa);
    setSolverAMGCL(prm, "solver.", sol_opt);

    /***************************************
     *        Solve problem                *
     ***************************************/
    auto t1 = std::chrono::high_resolution_clock::now();

    int block_size = 1;
    const int bz = 3;

    const auto matrix = amgcl::adapter::zero_copy(n, &cols[0], &rows[0], &entries[0]);
    switch(block_size){
      case 1:
      {
        amgcl::make_solver<
            amgcl::runtime::preconditioner<Backend>,
            amgcl::runtime::solver::wrapper<Backend>
        > solve(* matrix, prm);

        auto t2 = std::chrono::high_resolution_clock::now();
        if(verbose){
            std::cout << "Solver setup took "
                      << (double)std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count()/1000.0
                      << " seconds\n";
        }
        std::tie(iters, error) = solve(b, x);

        if(verbose){
            std::cout << solve << std::endl;
        }
      }
      case bz:
      {
        Backend::params bprm;
        typedef amgcl::backend::builtin<amgcl::static_matrix<double, bz, bz> > BBackend;
        amgcl::make_block_solver<
            amgcl::runtime::preconditioner<BBackend>,
            amgcl::runtime::solver::wrapper<BBackend>
        > solve(*matrix, prm, bprm);


        // auto t2 = std::chrono::high_resolution_clock::now();
        // if(verbose){
        //     std::cout << "Solver setup took "
        //               << (double)std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count()/1000.0
        //               << " seconds\n";
        // }
        // std::tie(iters, error) = solve(b, x);
        //
        // if(verbose){
        //     std::cout << solve << std::endl;
        // }
      }
    }
}

/* MEX gateway */

void mexFunction( int nlhs, mxArray *plhs[],
		  int nrhs, const mxArray *prhs[] )

{
    double *result;
    double *rhs;
    double *err;
    double *it_count;
    mwSize m,n,nnz;
    mwIndex * cols;
    mwIndex * rows;
    const mxArray * pa;

    double * entries;
    std::string relaxParam;
    std::string coarsenParam;

    if (nrhs != 6) {
	    mexErrMsgTxt("6 input arguments required.");
    } else if (nlhs > 3) {
	    mexErrMsgTxt("Wrong number of output arguments.");
    }

    m = mxGetM(prhs[0]);
    n = mxGetN(prhs[0]);
    int M = (int)m;

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
    plhs[2] = mxCreateDoubleMatrix(1, 1, mxREAL);

    result = mxGetPr(plhs[0]);
    err = mxGetPr(plhs[1]);
    it_count = mxGetPr(plhs[2]);

    cols    = mxGetJc(prhs[0]);
    rows    = mxGetIr(prhs[0]);
    entries = mxGetPr(prhs[0]);
    nnz  = mxGetNzmax(prhs[0]);
    rhs     = mxGetPr(prhs[1]);
    pa = prhs[2];
    double tolerance = mxGetScalar(prhs[3]);
    int maxiter = mxGetScalar(prhs[4]);
    int solver_strategy_id = mxGetScalar(prhs[5]);


    std::vector<double> b(n);
    for(int ix = 0; ix < n; ix++){
        b[ix] = rhs[ix];
    }
    int    iters;
    double error;
    std::vector<double> x(M, 0.0);
    switch(solver_strategy_id) {
        case 1:
            solve_regular(M, cols, rows, entries, pa, b, x, tolerance, maxiter, iters, error);
            break;
        case 2:
            solve_cpr(M, cols, rows, entries, pa, b, x, tolerance, maxiter, iters, error);
            break;
        default : mexErrMsgTxt("Unknown solver_strategy_id.");
    }

    for(int ix=0; ix < M; ix++){
        result[ix] = x[ix];
    }
    x.clear();
    b.clear();
    err[0] = error;
    it_count[0] = iters;
    return;
}
