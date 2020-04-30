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
// #include "solve_template.cpp"

/* Block system support */
#ifndef AMGCL_BLOCK_SIZES
#  define AMGCL_BLOCK_SIZES (2)(3)(4)(5)(6)(7)(8)(9)(10)
#endif
// #include "amgcl_block_macros.cpp"

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
typedef mwIndex indexType;

#define SOLVEBLOCK(bz) solve_block_system<bz>(n, tolerance, maxiter, mex_options, solver_type, ptr, col, V, rhs, result)

template <int B, typename matrix_type, typename rhs_type>
std::tuple<size_t, double> solve_regular(double tolerance, 
                                         int max_iter, 
                                         const mxArray* mex_options,
                                         const matrix_type matrix, 
                                         std::vector<rhs_type> &rhs, 
                                         std::vector<rhs_type> &x) {
    typedef amgcl::static_matrix<double, B, B> val_type;
    typedef amgcl::backend::builtin<val_type> Backend;

    // Get parameters from struct
    int relax_id = mxGetScalar(mxGetField(mex_options, 0, "relaxation"));
    bool verbose = mxGetScalar(mxGetField(mex_options, 0, "verbose"));
    int precond_id = mxGetScalar(mxGetField(mex_options, 0, "preconditioner"));

    std::string relaxParam;
    /***************************************
     *   Build parameter tree for solver   *
     ***************************************/

    boost::property_tree::ptree prm;

    // Set tolerance, max iterations and select preconditioner style
    prm.put("solver.tol", tolerance);
    if (max_iter > 0) {
        prm.put("solver.maxiter", max_iter);
    }
    switch (precond_id) {
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
    default: mexErrMsgTxt("Unknown precond_id.");
    }

    if (precond_id == 1) {
        // Select coarsening strategy
        amg_opts c_opt;
        setCoarseningStructMex(c_opt, mex_options);
        setCoarseningAMGCL(prm, "precond.", c_opt);
    }
    // Select relaxation strategy for solver
    relax_opts pr_opt;
    // pr_opt.relax_id = relax_id;
    setRelaxationStructMex(pr_opt, mex_options, "");
    setRelaxationAMGCL(prm, relaxParam, pr_opt);

    // Select solver
    solver_opts sol_opt;
    setSolverStructMex(sol_opt, mex_options);
    setSolverAMGCL(prm, "solver.", sol_opt);

    if (verbose) {
        std::cout << "Solver: Standard AMGCL solver." << std::endl;
    }
    // Typedef solver
    typedef amgcl::make_solver<
        amgcl::runtime::preconditioner<Backend>,
        amgcl::runtime::solver::wrapper<Backend>
    > Solver;
    // Instansiate and solve!
    Solver solve(matrix, prm); // Setup
    std::tuple<size_t, double> out;
    out = solve(rhs, x); // Solve
    if (verbose) {
        std::cout << solve << std::endl;
    }
    return out;
}

template <int B, typename matrix_type, typename rhs_type>
std::tuple<size_t, double> solve_cpr(double tolerance,
    int max_iter,
    const mxArray* mex_options,
    const matrix_type matrix,
    std::vector<rhs_type>& rhs,
    std::vector<rhs_type>& x) {
    typedef amgcl::static_matrix<double, B, B> val_type;
    typedef amgcl::backend::builtin<val_type> Backend;
    typedef amgcl::backend::builtin<double> ScalarBackend;
    boost::property_tree::ptree prm;
    // CPR settings
    int active_rows = mxGetScalar(mxGetField(mex_options, 0, "active_rows"));
    bool use_drs = mxGetScalar(mxGetField(mex_options, 0, "use_drs"));
    // Pressure and global relaxation choices
    int relax_p_id = mxGetScalar(mxGetField(mex_options, 0, "relaxation"));
    int relax_s_id = mxGetScalar(mxGetField(mex_options, 0, "s_relaxation"));
    // Various settings
    bool verbose = mxGetScalar(mxGetField(mex_options, 0, "verbose"));
    /*****************************************
     * Begin building parameter tree for CPR *
     ****************************************/
    /* Set tolerance */
    prm.put("solver.tol", tolerance);
    if (max_iter > 0) {
        prm.put("solver.maxiter", max_iter);
    }
    prm.put("precond.block_size", B);
    prm.put("precond.active_rows", active_rows);
    /* Select coarsening strategy */
    amg_opts c_opt;
    setCoarseningStructMex(c_opt, mex_options);
    setCoarseningAMGCL(prm, "precond.pprecond.", c_opt);

    /* Select relaxation strategy for pressure solver */
    relax_opts pr_opt;
    // pr_opt.relax_id = relax_p_id;
    setRelaxationStructMex(pr_opt, mex_options, "");
    setRelaxationAMGCL(prm, "precond.pprecond.relax.", pr_opt);
    /* Select relaxation strategy for second stage solver */
    relax_opts ps_opt;
    // ps_opt.relax_id = relax_s_id;
    setRelaxationStructMex(ps_opt, mex_options, "s_");
    setRelaxationAMGCL(prm, "precond.sprecond.", ps_opt);

    /* Select solver */
    solver_opts sol_opt;
    setSolverStructMex(sol_opt, mex_options);
    setSolverAMGCL(prm, "solver.", sol_opt);

    // Output iterations and residual norm
    std::tuple<size_t, double> out;

    // Pressure solver for CPR
    typedef amgcl::amg<ScalarBackend,
        amgcl::runtime::coarsening::wrapper,
        amgcl::runtime::relaxation::wrapper
    > PPrecond;
    // Second-stage solver for CPR
    typedef amgcl::relaxation::as_preconditioner<Backend, amgcl::runtime::relaxation::wrapper>
        SPrecond;
    if (use_drs) {
        if (verbose) {
            std::cout << "Solver: Dynamic-row-sum CPR" << std::endl;
        }
        double dd = mxGetScalar(mxGetField(mex_options, 0, "drs_eps_dd"));
        double ps = mxGetScalar(mxGetField(mex_options, 0, "drs_eps_ps"));
        prm.put("precond.eps_dd", dd);
        prm.put("precond.eps_ps", ps);

        mxArray* drs_weights_mx = mxGetField(mex_options, 0, "drs_row_weights");
        size_t drs_weights_n = mxGetM(drs_weights_mx);

        if (drs_weights_n > 0) {
            double* drs_weights = mxGetPr(drs_weights_mx);
            prm.put("precond.weights", drs_weights);
            prm.put("precond.weights_size", drs_weights_n);
        }
        typedef amgcl::make_solver<
            amgcl::preconditioner::cpr_drs<PPrecond, SPrecond>,
            amgcl::runtime::solver::wrapper<Backend>
        > Solver;

        // Instansiate and solve!
        Solver solve(matrix, prm); // Setup

        out = solve(rhs, x); // Solve

        if (verbose) {
            std::cout << solve << std::endl;
        }
    }else{
        if (verbose) {
            std::cout << "Solver: Regular CPR" << std::endl;
        }
        // Regular CPR
        typedef amgcl::make_solver<
            amgcl::preconditioner::cpr<PPrecond, SPrecond>,
            amgcl::runtime::solver::wrapper<Backend>
        > Solver;

        // Instansiate and solve!
        Solver solve(matrix, prm); // Setup
        
        out = solve(rhs, x); // Solve

        if (verbose) {
            std::cout << solve << std::endl;
        }
    }
    return out;
}

template <int B>
std::tuple<size_t, double> solve_block_system(int n, double tolerance, int max_iter, const mxArray* mex_options, const int solver_id, indexType * ptr, indexType * col, double * V, double * rhs, double * result) {
    typedef amgcl::static_matrix<double, B, B> val_type;
    typedef amgcl::static_matrix<double, B, 1> rhs_type;

    const val_type* v_ptr = reinterpret_cast<const val_type*>(V);
    const auto matrix = amgcl::adapter::zero_copy(n, ptr, col, v_ptr);

    /* We have the matrix and solver, deal with input and right-hand side */

    std::vector<rhs_type> b(n);
#pragma omp parallel for schedule(static)
    for (int i = 0; i < n; i++) {
        auto tmp = b[i].data();
        for (int j = 0; j < B; j++) {
            // mexPrintf("El %d block %d: rhs %f\n", i, j, rhs[i * block_size + j]);
            // tmp[j] = rhs[i * B + j];
            tmp[j] = rhs[i + j*n];
        }
    }

    std::vector<rhs_type> x(n);

    // std::cout << "Tolerance: " << tolerance << " Max iter: " << max_iter << std::endl;
    /* Perform the solve */
    std::tuple<size_t, double> out;
    switch (solver_id) {
        case 1:
            // Regular solver
            out = solve_regular<B>(tolerance, max_iter, mex_options, matrix, b, x);
            break;
        case 2:
            out = solve_cpr<B>(tolerance, max_iter, mex_options, matrix, b, x);
            break;
        default:
            mexErrMsgTxt("Block size not supported!");
    }
    
    // mexPrintf("Used %d iterations for error of %g\n", iters, error);

    /* Copy the std::vector block buffer into the scalar MEX array C-style raw buffer */
#pragma omp parallel for schedule(static)
    for (int i = 0; i < n; i++) {
        auto tmp = x[i];
        for (int j = 0; j < B; j++) {
            // result[i * B + j] = tmp(j);
            result[i + j*n] = tmp(j);
        }
    }
    return out;
}


/* MEX gateway */
void mexFunction(int nlhs, mxArray* plhs[],
    int nrhs, const mxArray* prhs[])

{
    double* result;
    double* err;
    double* it_count;
    mwSize m, n, nnz, m_rhs, n_rhs;
    mwIndex* cols;
    mwIndex* rows;
    const mxArray* pa;

    double* entries;
    std::string relaxParam;
    std::string coarsenParam;

    if (nrhs == 0 && nlhs == 0) {
        mexPrintf("AMGCL is compiled and ready for use.\n");
        return;
    }
    else if (nrhs != 7 && nrhs != 8) {
        mexErrMsgTxt("7 or 8 input arguments required.\nSyntax: amgcl_matlab(I, J, V, b, amg_opt, tol, maxit, (id optional))");
    }
    else if (nlhs > 3) {
        mexErrMsgTxt("More than three outputs requested!");
    }
    /*
    *  Input argument & validation.
    */

    // First three inputs: Matrix pointers
    n = mxGetNumberOfElements(prhs[1]) - 1; // rows/cols in matrix
    size_t block_prod = mxGetM(prhs[2]); // block size
    size_t block_size = sqrt(block_prod);
    int M = n * block_size;
    // OK!
    
    indexType* I = (indexType*)mxGetData(prhs[0]);
    indexType* J = (indexType*)mxGetData(prhs[1]);

    indexType* col = I;
    indexType* ptr = J;

    int solver_type;
    if (nrhs < 8) {
        solver_type = 1;
    }
    else {
        solver_type = mxGetScalar(prhs[7]);
    }

    // Third argument: Matrix of block-entries with each block sequentially in memory
    double * V = mxGetPr(prhs[2]);
    int n_elements = mxGetN(prhs[2]);
    // Fourth argument: Right-hand side.
    double * rhs = mxGetPr(prhs[3]);

    // Five: MEX options struct
    const mxArray* mex_options = prhs[4];

    // Six and seven: Tolerance and maximum number of iterations
    double tolerance = mxGetScalar(prhs[5]);
    int maxiter = (int)mxGetScalar(prhs[6]);

    // mexPrintf("Recieved a matrix of size %d with %d elements and block size %d x %d\n", n, n_elements, block_size, block_size);


    /* 
    *  Output arguments.
    */

    // First output: Solution column vector
    plhs[0] = mxCreateDoubleMatrix(M, 1, mxREAL);
    // Second output: Residual error
    plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
    // Third output: Number of iterations
    plhs[2] = mxCreateDoubleMatrix(1, 1, mxREAL);

    // Get pointers to output
    result = mxGetPr(plhs[0]);
    err = mxGetPr(plhs[1]);
    it_count = mxGetPr(plhs[2]);

    int    iters = 0;
    double error = 0;
    const size_t B = 2;

    /*
    *  Solving
    */

    if (0) {
        int outer = 0;
        for (int rowNo = 0; rowNo < n; rowNo++) {
            for (int pos = J[rowNo]; pos < J[rowNo + 1]; pos++) {
                int colNo = I[pos];
                mexPrintf("Element %d:\n* %d, %d:\n", outer+1, rowNo+1, colNo+1);
                for (int bi = 0; bi < block_size; bi++) {
                    for (int bj = 0; bj < block_size; bj++) {
                        mexPrintf("%5g ", V[outer*block_prod + bi*block_size + bj]);
                    }
                    mexPrintf("\n");
                }
                outer++;
            }
        }
    }
    switch(block_size){
        case 2:
            std::tie(iters, error) = SOLVEBLOCK(2);
            break;
        case 3:
            std::tie(iters, error) = SOLVEBLOCK(3);
            break;
        case 4:
            std::tie(iters, error) = SOLVEBLOCK(4);
            break;
        case 5:
            std::tie(iters, error) = SOLVEBLOCK(5);
            break;
        case 6:
            std::tie(iters, error) = SOLVEBLOCK(6);
            break;
        case 7:
            std::tie(iters, error) = SOLVEBLOCK(7);
            break;
        case 8:
            std::tie(iters, error) = SOLVEBLOCK(8);
            break;
        case 9:
            std::tie(iters, error) = SOLVEBLOCK(9);
            break;
        case 10:
            std::tie(iters, error) = SOLVEBLOCK(10);
            break;
        default:
            mexErrMsgTxt("Block size not supported!");
    }
    err[0] = error;
    it_count[0] = iters;
    return;
}
