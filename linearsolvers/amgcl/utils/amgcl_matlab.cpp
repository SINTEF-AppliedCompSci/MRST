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



/* Relaxation */
struct relax_opts {
    int relax_id;
};

void setRelaxationAMGCL(boost::property_tree::ptree & prm, std::string relaxParam, relax_opts opts){
    std::string relaxType = relaxParam + "type";
    switch(opts.relax_id) {
        case 1: 
            prm.put(relaxType,  amgcl::runtime::relaxation::spai0);
            break;
        case 2: 
            prm.put(relaxType,  amgcl::runtime::relaxation::gauss_seidel);
            break;
        case 3: 
            prm.put(relaxType,  amgcl::runtime::relaxation::ilu0);
            break;
        case 4: 
            prm.put(relaxType,  amgcl::runtime::relaxation::iluk);
            break;
        case 5: 
            prm.put(relaxType,  amgcl::runtime::relaxation::ilut);
            break;
        case 6: 
            prm.put(relaxType,  amgcl::runtime::relaxation::damped_jacobi);
            break;
        case 7: 
            prm.put(relaxType,  amgcl::runtime::relaxation::spai1);
            break;
        case 8: 
            prm.put(relaxType,  amgcl::runtime::relaxation::chebyshev);
            break;
        default : mexErrMsgTxt("Unknown relax_id."); 
    }
}

/* Coarsening */ 
struct amg_opts {
    int coarsen_id;
    int coarse_enough;
    bool direct_coarse;
    int max_levels;
    int ncycle;
    int npre;
    int npost;
    int pre_cycles;
};

void setCoarseningStructMex(amg_opts &c_opt, const mxArray * pa){
    /* Convert mex struct pointer to struct for amg */
    c_opt.coarsen_id = mxGetScalar(mxGetField(pa, 0, "coarsening"));
    c_opt.coarse_enough = mxGetScalar(mxGetField(pa, 0, "coarse_enough"));
    c_opt.direct_coarse = mxGetScalar(mxGetField(pa, 0, "direct_coarse"));
    c_opt.max_levels = mxGetScalar(mxGetField(pa, 0, "max_levels"));
    c_opt.ncycle = mxGetScalar(mxGetField(pa, 0, "ncycle"));
    c_opt.npre = mxGetScalar(mxGetField(pa, 0, "npre"));
    c_opt.npost = mxGetScalar(mxGetField(pa, 0, "npost"));
    c_opt.pre_cycles = mxGetScalar(mxGetField(pa, 0, "pre_cycles"));
}

void setCoarseningAMGCL(boost::property_tree::ptree & prm, std::string prefix, amg_opts options){
    std::string coarsetype = prefix + "coarsening.type";
    switch(options.coarsen_id) {
        case 1: 
            prm.put(coarsetype,  amgcl::runtime::coarsening::smoothed_aggregation);
            break;
        case 2: 
            prm.put(coarsetype,  amgcl::runtime::coarsening::ruge_stuben);
            break;
        case 3: 
            prm.put(coarsetype,  amgcl::runtime::coarsening::aggregation);
            break;
        case 4: 
            prm.put(coarsetype,  amgcl::runtime::coarsening::smoothed_aggr_emin);
            break;
        default : mexErrMsgTxt("Unknown coarsen_id: " + options.coarsen_id); 
    }
    /* When is a level coarse enough */
    if (options.coarse_enough >= 0){
        prm.put(prefix + "coarse_enough", options.coarse_enough);
    }
    /* Use direct solver for coarse sys */
    prm.put(prefix + "direct_coarse", options.direct_coarse);
    /* Max levels */
    if (options.max_levels >= 0){
        prm.put(prefix + "max_levels", options.max_levels);
    }
    /* Number of cycles */
    if (options.ncycle >= 0){
        prm.put(prefix + "ncycle", options.ncycle);
    }
    /* Pre cycles */
    if (options.npre >= 0){
        prm.put(prefix + "npre", options.npre);
    }
    /* Post cycles */
    if (options.npost >= 0){
        prm.put(prefix + "npost", options.npost);
    }
    /* Pre cycles (precond) */
    if (options.pre_cycles >= 0){
        prm.put(prefix + "pre_cycles", options.pre_cycles);
    }
}
/* Krylov solver */
struct solver_opts {
    int solver_id;
    int L;
    int M;
    int K;
    int S;
    double delta;
    double omega;
    bool convex;
    bool always_reset;
    bool store_Av;
    bool replace;
};

void setSolverStructMex(solver_opts &opt, const mxArray * pa){
    opt.solver_id = mxGetScalar(mxGetField(pa, 0, "solver"));
    opt.L = mxGetScalar(mxGetField(pa, 0, "bicgstabl_l"));
    opt.M = mxGetScalar(mxGetField(pa, 0, "gmres_m"));
    opt.K = mxGetScalar(mxGetField(pa, 0, "lgmres_k"));
    opt.S = mxGetScalar(mxGetField(pa, 0, "idrs_s"));
    opt.delta = mxGetScalar(mxGetField(pa, 0, "bicgstabl_delta"));
    opt.omega = mxGetScalar(mxGetField(pa, 0, "idrs_omega"));
    opt.convex = mxGetScalar(mxGetField(pa, 0, "bicgstabl_convex"));
    opt.always_reset = mxGetScalar(mxGetField(pa, 0, "lgmres_always_reset"));
    opt.store_Av = mxGetScalar(mxGetField(pa, 0, "lgmres_store_av"));
    opt.replace = mxGetScalar(mxGetField(pa, 0, "idrs_replacement"));
}


void setSolverAMGCL(boost::property_tree::ptree & prm, std::string prefix, solver_opts options){
    std::string solvertype = prefix + "type";
    switch(options.solver_id) {
        case 1: 
            prm.put(solvertype,  amgcl::runtime::solver::bicgstab);
            break;
        case 2: 
            prm.put(solvertype,  amgcl::runtime::solver::cg);
            break;
        case 3: 
            prm.put(solvertype,  amgcl::runtime::solver::bicgstabl);
            {
                prm.put(prefix + "L", options.L);
                prm.put(prefix + "delta", options.delta);
                prm.put(prefix + "convex", options.convex);
            }
            break;
        case 4: 
            prm.put(solvertype,  amgcl::runtime::solver::gmres);
            {
                prm.put(prefix + "M", options.M);
            }
            break;
        case 5: 
            prm.put(solvertype,  amgcl::runtime::solver::lgmres);
            {
                prm.put(prefix + "M", options.M);
                prm.put(prefix + "K", options.K);
                prm.put(prefix + "always_reset", options.always_reset);
                prm.put(prefix + "store_Av", options.store_Av);            
            }
            break;
        case 6: 
            prm.put(solvertype,  amgcl::runtime::solver::fgmres);
            {
                prm.put(prefix + "M", options.M);
            }
            break;
        case 7: 
            prm.put(solvertype,  amgcl::runtime::solver::idrs);
            {
                prm.put(prefix + "s", options.S);
                prm.put(prefix + "omega", options.omega);
                prm.put(prefix + "replacement", options.replace);
            }
            break;
        default : mexErrMsgTxt("Unknown solver_id."); 
    }
}


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
    prm.put("precond.block_size", block_size);
    prm.put("precond.active_rows", active_rows);
    /* Select coarsening strategy */
    amg_opts c_opt;
    setCoarseningStructMex(c_opt, pa);
    setCoarseningAMGCL(prm, "precond.pprecond.", c_opt);

    /* Select relaxation strategy for pressure solver */
    relax_opts pr_opt;
    pr_opt.relax_id = relax_p_id;
    setRelaxationAMGCL(prm, "precond.pprecond.relax.", pr_opt);
    /* Select relaxation strategy for second stage solver */
    relax_opts ps_opt;
    ps_opt.relax_id = relax_s_id;
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
    pr_opt.relax_id = relax_id;
    setRelaxationAMGCL(prm, relaxParam, pr_opt);
    
    /* Select solver */
    solver_opts sol_opt;
    setSolverStructMex(sol_opt, pa);
    setSolverAMGCL(prm, "solver.", sol_opt);

    /***************************************
     *        Solve problem                *
     ***************************************/
    auto t1 = std::chrono::high_resolution_clock::now();


    amgcl::make_solver<
        amgcl::runtime::preconditioner<Backend>,
        amgcl::runtime::solver::wrapper<Backend>
    > solve(amgcl::adapter::zero_copy(n, &cols[0], &rows[0], &entries[0]), prm);
    
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


