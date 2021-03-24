//
// include necessary system headers
//
#define _USE_MATH_DEFINES

#include <cmath>
#include <mex.h>
#include <array>
#include <memory>
#include <string>
#include <chrono>
#include <iostream>

#ifndef HAVE_OCTAVE
#include "matrix.h"
#endif

#include <amgcl/make_solver.hpp>
#include <amgcl/backend/builtin.hpp>
// AMG, relaxation etc
#include <amgcl/amg.hpp>
#include <amgcl/relaxation/as_preconditioner.hpp>
// Matrix adapters
#include <amgcl/adapter/crs_tuple.hpp>
#include <amgcl/backend/block_crs.hpp>
#include <amgcl/adapter/zero_copy.hpp>
// Include runtime parameters
#include <amgcl/coarsening/runtime.hpp>
#include <amgcl/relaxation/runtime.hpp>
#include <amgcl/solver/runtime.hpp>
#include <amgcl/preconditioner/runtime.hpp>

// CPR
#include <amgcl/preconditioner/cpr.hpp>
#include <amgcl/preconditioner/cpr_drs.hpp>
// Utilities etc
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/preprocessor/seq/for_each.hpp>
#include <boost/preprocessor/cat.hpp>

/* MEX interfaces */
#include "amgcl_mex_utils.cpp"
#include "solve_template.cpp"

/* Block system support */
#ifndef AMGCL_BLOCK_SIZES
#  define AMGCL_BLOCK_SIZES (2)(3)(4)(5)(6)(7)(8)(9)(10)
#endif
#include "amgcl_block_macros.cpp"

#ifndef SOLVER_BACKEND_BUILTIN
#  define SOLVER_BACKEND_BUILTIN
#endif

#include <amgcl/backend/builtin.hpp>
#include <amgcl/value_type/static_matrix.hpp>
#include <amgcl/adapter/block_matrix.hpp>
#include <amgcl/make_block_solver.hpp>
#include <amgcl/adapter/crs_tuple.hpp>

#ifdef _OPENMP
   #include <omp.h>
#else
   #define omp_get_max_threads() 1
#endif

typedef amgcl::backend::builtin<double> Backend;
// Scalar solver
typedef amgcl::make_solver<
    amgcl::runtime::preconditioner<Backend>,
    amgcl::runtime::solver::wrapper<Backend>
> ScalarSolver;

// Shared pointer for scalar solver
static std::shared_ptr<ScalarSolver> scalar_solve_ptr(nullptr);

// Pressure solver for CPR
typedef amgcl::amg<Backend,
    amgcl::runtime::coarsening::wrapper,
    amgcl::runtime::relaxation::wrapper
> PPrecond;
// Second-stage solver for CPR
typedef amgcl::relaxation::as_preconditioner<Backend, amgcl::runtime::relaxation::wrapper>
    SPrecond;
// Regular CPR
typedef amgcl::make_solver<
            amgcl::preconditioner::cpr<PPrecond, SPrecond>,
            amgcl::runtime::solver::wrapper<Backend>
            > CPRSolver;

static std::shared_ptr<CPRSolver> cpr_solve_ptr(nullptr);

// CPR with dynamic row sum
typedef amgcl::make_solver<
            amgcl::preconditioner::cpr_drs<PPrecond, SPrecond>,
            amgcl::runtime::solver::wrapper<Backend>
            > CPRSolverDRS;

static std::shared_ptr<CPRSolverDRS> cpr_drs_solve_ptr(nullptr);

BOOST_PP_SEQ_FOR_EACH(AMGCL_DEFINE_BLOCK_TYPES, ~, AMGCL_BLOCK_SIZES)
BOOST_PP_SEQ_FOR_EACH(AMGCL_DEFINE_BLOCK_SOLVER, BlockSolverSize, AMGCL_BLOCK_SIZES)
BOOST_PP_SEQ_FOR_EACH(AMGCL_DEFINE_BLOCK_CPR_SOLVERS, ~, AMGCL_BLOCK_SIZES)

static void reset_solvers(void){
    // Reset scalar solver
    scalar_solve_ptr.reset();
    // Reset CPR
    cpr_solve_ptr.reset();
    // Reset CPR-DRS
    cpr_drs_solve_ptr.reset();
    // Reset basic block solvers
    BOOST_PP_SEQ_FOR_EACH(AMGCL_RESET_BLOCK_SOLVER, block_solve_ptr, AMGCL_BLOCK_SIZES)
    // Reset CPR variants
    BOOST_PP_SEQ_FOR_EACH(AMGCL_RESET_BLOCK_SOLVER, cpr_block_solve_ptr, AMGCL_BLOCK_SIZES)
    BOOST_PP_SEQ_FOR_EACH(AMGCL_RESET_BLOCK_SOLVER, cpr_drs_block_solve_ptr, AMGCL_BLOCK_SIZES)
}
// CPR Gateway
template <class M>
void solve_cpr(int n, const M matrix, const mxArray * pa,
        std::vector<double> b, std::vector<double> & x, double tolerance,
        int maxiter, int & iters, double & error){

    // CPR settings
    // bool update_s   = mxGetScalar(mxGetField(pa, 0, "update_sprecond"));
    bool update_s = GET_STRUCT_SCALAR(pa, "update_sprecond");
    bool update_p   = GET_STRUCT_SCALAR(pa, "update_ptransfer");
    bool use_blocks = GET_STRUCT_SCALAR(pa, "cpr_blocksolver");
    int block_size  = GET_STRUCT_SCALAR(pa, "block_size");
    int active_rows = GET_STRUCT_SCALAR(pa, "active_rows");
    bool use_drs    = GET_STRUCT_SCALAR(pa, "use_drs");
    
    // Pressure and global relaxation choices
    int relax_p_id  = GET_STRUCT_SCALAR(pa, "relaxation");
    int relax_s_id  = GET_STRUCT_SCALAR(pa, "s_relaxation");
    // Various settings
    bool verbose      = GET_STRUCT_SCALAR(pa, "verbose");
    bool write_params = GET_STRUCT_SCALAR(pa, "write_params");
    /*****************************************
     * Begin building parameter tree for CPR *
     ****************************************/
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
    if(use_drs){
        double dd = GET_STRUCT_SCALAR(pa, "drs_eps_dd");
        double ps = GET_STRUCT_SCALAR(pa, "drs_eps_ps");
        prm.put("precond.eps_dd", dd);
        prm.put("precond.eps_ps", ps);

        mxArray * drs_weights_mx = mxGetField(pa, 0, "drs_row_weights");
        size_t drs_weights_n = mxGetM(drs_weights_mx);

        if(drs_weights_n>0){
            double * drs_weights = mxGetPr(drs_weights_mx);
            prm.put("precond.weights", drs_weights);
            prm.put("precond.weights_size", drs_weights_n);
        }
        if(write_params){
            std::cout << "Writing amgcl setup file to mrst_amgcl_cpr_drs_setup.json" << std::endl;
            std::ofstream file("mrst_amgcl_drs_setup.json");
            boost::property_tree::json_parser::write_json(file, prm);
        }
        if(!use_blocks){
          std::tie(iters, error) = solve_shared_cpr(cpr_drs_solve_ptr, *matrix, b, x, prm, matrix->nrows, update_s, update_p, verbose);
        }else{
          switch(block_size){
            BOOST_PP_SEQ_FOR_EACH(AMGCL_BLOCK_CPR_SOLVER, cpr_drs_block_solve_ptr, AMGCL_BLOCK_SIZES)
            default:
                mexErrMsgIdAndTxt("AMGCL:UndefBlockSize",
                                  "Failure: Block size %d not supported.",
                                  block_size);
          }
        }
    }else{
         if(write_params){
            std::cout << "Writing amgcl setup file to mrst_amgcl_cpr_setup.json" << std::endl;
            std::ofstream file("mrst_amgcl_setup.json");
            boost::property_tree::json_parser::write_json(file, prm);
        }
        if(!use_blocks){
          std::tie(iters, error) = solve_shared_cpr(cpr_solve_ptr, *matrix, b, x, prm, matrix->nrows, update_s, update_p, verbose);
        }else{
          switch(block_size){
            BOOST_PP_SEQ_FOR_EACH(AMGCL_BLOCK_CPR_SOLVER, cpr_block_solve_ptr, AMGCL_BLOCK_SIZES)
            default:
                mexErrMsgIdAndTxt("AMGCL:UndefBlockSize",
                                  "Failure: Block size %d not supported.",
                                  block_size);
          }
        }
    }
}

template <class M>
void solve_regular(int n, const M matrix, const mxArray * pa,
        const std::vector<double> & b, std::vector<double> & x, double tolerance,
        int maxiter, int & iters, double & error){
    // Get parameters from struct
    int relax_id      = GET_STRUCT_SCALAR(pa, "relaxation");
    bool verbose      = GET_STRUCT_SCALAR(pa, "verbose");
    bool write_params = GET_STRUCT_SCALAR(pa, "write_params");
    int precond_id    = GET_STRUCT_SCALAR(pa, "preconditioner");
    std::string relaxParam;
    /***************************************
     *   Build parameter tree for solver   *
     ***************************************/

    boost::property_tree::ptree prm;
    /* Set tolerance, max iterations and select preconditioner style */
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
    if(write_params){
      std::cout << "Writing amgcl setup file to mrst_amgcl_cpr_setup.json" << std::endl;
      std::ofstream file("mrst_regular_setup.json");
      boost::property_tree::json_parser::write_json(file, prm);
    }
    int block_size = GET_STRUCT_SCALAR(pa, "block_size");
    switch(block_size){
      case 0:
      case 1:
      {
        std::tie(iters, error) = solve_shared(scalar_solve_ptr, matrix, b, x, prm, verbose);
      } break;
      BOOST_PP_SEQ_FOR_EACH(AMGCL_BLOCK_SOLVER, block_solve_ptr, AMGCL_BLOCK_SIZES)
        default:
            mexErrMsgIdAndTxt("AMGCL:UndefBlockSize",
                              "Failure: Block size %d not supported.",
                              block_size);
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
    mwSize m, n, nnz, m_rhs, n_rhs, m_initial_guess, n_initial_guess;
    mwIndex * cols;
    mwIndex * rows;
    const mxArray * pa;

    double * entries;
    std::string relaxParam;
    std::string coarsenParam;
    
    if (nrhs == 0 && nlhs == 0){
        mexPrintf("AMGCL is compiled and ready for use.\n");
        return;
    } else if (nrhs != 6 && nrhs != 7 && nrhs != 8) {
	    mexErrMsgTxt("6, 7 or 8 input arguments required.\nSyntax: amgcl_matlab(A, b, opts, tol, maxit, solver_id, reuse_id, x0)");
    } else if (nlhs > 3) {
	    mexErrMsgTxt("More than three outputs requested!");
    }

    m = mxGetM(prhs[0]);
    n = mxGetN(prhs[0]);
    
    m_rhs = mxGetM(prhs[1]);
    n_rhs = mxGetN(prhs[1]);
    
    bool has_initial_guess;
    double *initial_guess;
    has_initial_guess = (nrhs == 8);
    if (has_initial_guess){
        m_initial_guess = mxGetM(prhs[7]);
        n_initial_guess = mxGetN(prhs[7]);
    }
    
    int M = (int)m;

    if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) ||  !mxIsSparse(prhs[0]) ) {
	    mexErrMsgTxt("Matrix should be a real sparse matrix.");
        return;
    }
    if (n != m) {
	    mexErrMsgTxt("Matrix must be square.");
        return;
    }
    if (!mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) || n_rhs != 1 || m_rhs != m) {
	    mexErrMsgTxt("Right hand side must be real double column vector with one entry per row of A.");
        return;
    }
    if (has_initial_guess && (!mxIsDouble(prhs[7]) || mxIsComplex(prhs[7]) || n_initial_guess != 1 || m_initial_guess != m)) {
	    mexErrMsgTxt("Initial guess must be real double column vector with one entry per column of A.");
        return;
    }
    // First output: Solution column vector
    plhs[0] = mxCreateDoubleMatrix(m, 1, mxREAL);
    // Second output: Residual error
    plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
    // Third output: Number of iterations
    plhs[2] = mxCreateDoubleMatrix(1, 1, mxREAL);

    result   = mxGetPr(plhs[0]);
    err      = mxGetPr(plhs[1]);
    it_count = mxGetPr(plhs[2]);

    cols          = mxGetJc(prhs[0]);
    rows          = mxGetIr(prhs[0]);
    entries       = mxGetPr(prhs[0]);
    nnz           = mxGetNzmax(prhs[0]);
    rhs           = mxGetPr(prhs[1]);
    if (has_initial_guess){
        initial_guess = mxGetPr(prhs[7]);
    }
    // Build system matrix
    const auto matrix = amgcl::adapter::zero_copy(n, &cols[0], &rows[0], &entries[0]);
    // Get struct with options
    pa = prhs[2];
    double tolerance       = mxGetScalar(prhs[3]);
    int maxiter            = (int)mxGetScalar(prhs[4]);
    int solver_strategy_id = (int)mxGetScalar(prhs[5]);
    bool verbose           = GET_STRUCT_SCALAR(pa, "verbose");
    int nthreads           = GET_STRUCT_SCALAR(pa, "nthreads");
    int block_size         = GET_STRUCT_SCALAR(pa, "block_size");
    int reuse_mode;
    if(nrhs == 7){
      reuse_mode = (int)mxGetScalar(prhs[6]);
    }else{
      reuse_mode = 1;
    }
    if(verbose){
      std::cout << "AMGCL solver recieved problem with " << n << " degrees of freedom.";
      if(block_size == 1){
        std::cout << " Treating system as scalar.";
      }else{
        std::cout << " System has block size of " << block_size << ".";
      }
      std::cout << std::endl;
      if(nthreads == 1){
        std::cout << "Solving in serial.";
      }
      else{
        std::cout << "Solving with " << nthreads << " threads using OpenMP.";
      }      
      std::cout << std::endl;
    }
    #ifdef _OPENMP
        omp_set_num_threads(nthreads);
    #endif
    std::vector<double> b(n);
    #pragma omp parallel for
    for(int ix = 0; ix < n; ix++){
        b[ix] = rhs[ix];
    }
    int    iters;
    double error;
    std::vector<double> x(M, 0.0);
    if (has_initial_guess){
        #pragma omp parallel for
        for(int ix = 0; ix < M; ix++){
            x[ix] = initial_guess[ix];
        }
    }
    switch(reuse_mode){
      case 1:
          // Default: No reuse, delete all if present
          reset_solvers();
          break;
      case 2:
          // Perform reuse
          break;
      default : mexErrMsgTxt("Unknown reuse mode: Must be 1 for no reuse or 2 for reuse.");
    }
    switch(solver_strategy_id) {
        case 1:
            solve_regular(M, matrix, pa, b, x, tolerance, maxiter, iters, error);
            break;
        case 2:
            solve_cpr(M, matrix, pa, b, x, tolerance, maxiter, iters, error);
            break;
        case 1000:
            // Remove shared pointers
            if(verbose){
                std::cout << "Resetting all solvers." << std::endl;
            }
            reset_solvers();
            break;
        default : mexErrMsgTxt("Unknown solver_strategy_id.");
    }
    if(reuse_mode == 1){
      reset_solvers();
    }
    #pragma omp parallel for
    for(int ix=0; ix < M; ix++){
        result[ix] = x[ix];
    }
    x.clear();
    b.clear();
    err[0] = error;
    it_count[0] = iters;
    mexAtExit(reset_solvers);
    return;
}
