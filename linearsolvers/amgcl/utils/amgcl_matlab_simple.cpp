//
// include necessary system headers
//
#define _USE_MATH_DEFINES
#include "matrix.h"

#include <cmath>
#include <mex.h>
#include <array>
#ifndef HAVE_OCTAVE
#include "matrix.h"
#endif
#include <memory>
#include <string>
#include <chrono>
#include <iostream>

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
//#include "amgcl_mex_utils.cpp"
#include "solve_template.cpp"

/* Block system support */
//(4)(5)(6)(7)(8)(9)(10)
#ifndef AMGCL_BLOCK_SIZES
#  define AMGCL_BLOCK_SIZES (2)(3)
//#  define AMGCL_BLOCK_SIZES (2)(3)//(4)(5)(6)(7)(8)(9)(10)
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
    // // Reset basic block solvers
    BOOST_PP_SEQ_FOR_EACH(AMGCL_RESET_BLOCK_SOLVER, block_solve_ptr, AMGCL_BLOCK_SIZES)
    // // Reset CPR variants
    BOOST_PP_SEQ_FOR_EACH(AMGCL_RESET_BLOCK_SOLVER, cpr_block_solve_ptr, AMGCL_BLOCK_SIZES)
    BOOST_PP_SEQ_FOR_EACH(AMGCL_RESET_BLOCK_SOLVER, cpr_drs_block_solve_ptr, AMGCL_BLOCK_SIZES)
}



void solve_cpr(int n, mwIndex * cols, mwIndex * rows, double * entries,
	       std::vector<double> b, std::vector<double> & x, double tolerance,
	       int maxiter, int & iters, double & error, boost::property_tree::ptree prm){
    /***************************************
     * Start AMGCL-link and select options *
     ***************************************/

    typedef
    amgcl::amg<Backend, amgcl::runtime::coarsening::wrapper, amgcl::runtime::relaxation::wrapper>
        PPrecond;

    typedef
    amgcl::relaxation::as_preconditioner<Backend, amgcl::runtime::relaxation::wrapper>
            SPrecond;
    //boost::property_tree::ptree prm;
    /* Set tolerance */
      
    /***************************************
     *        Solve problem                *
     ***************************************/
    auto t1 = std::chrono::high_resolution_clock::now();
    const auto matrix = amgcl::adapter::zero_copy(n, &cols[0], &rows[0], &entries[0]);
    bool use_drs = prm.get<bool>("use_drs");
    bool use_blocks = prm.get<bool>("cpr_blocksolver");
    bool verbose = prm.get<int>("verbosity")>0;
    bool write_params = prm.get<bool>("write_params");
    bool update_s   =  prm.get<bool>("update_sprecond");
    bool update_p  =  prm.get<bool>("update_ptransfer");
    int  block_size = prm.get<int>("block_size");
    if(use_drs){   
        if(prm.get<int>("verbosity")>10){
            std::cout << "Writing amgcl setup file to mrst_amgcl_cpr_drs_setup.json" << std::endl;
            std::ofstream file("mrst_amgcl_drs_setup.json");
            boost::property_tree::json_parser::write_json(file, prm);
        }
        
        if(!use_blocks){
          std::tie(iters, error) = solve_shared_cpr(cpr_drs_solve_ptr, *matrix, b, x, prm, matrix->nrows, update_s, update_p, verbose);
          if(verbose){
              std::cout << *cpr_drs_solve_ptr << std::endl;
          }
        }else{
          switch(block_size){
            BOOST_PP_SEQ_FOR_EACH(AMGCL_BLOCK_CPR_SOLVER, cpr_drs_block_solve_ptr, AMGCL_BLOCK_SIZES)
            default:
                mexErrMsgIdAndTxt("AMGCL:UndefBlockSize",
                                  "Failure: Block size %d not supported.",
                                  block_size);
          }
          // if(verbose){
          //     std::cout << cpr_drs_block_solve_ptr << std::endl;
          // }
        }

    }else{
         bool write_params = prm.get<bool>("write_params");
         if(write_params){
            std::cout << "Writing amgcl setup file to mrst_amgcl_cpr_setup.json" << std::endl;
            std::ofstream file("mrst_amgcl_setup.json");
            boost::property_tree::json_parser::write_json(file, prm);
        }
        
	if(!use_blocks){
	    std::tie(iters, error) = solve_shared_cpr(cpr_solve_ptr, *matrix, b, x, prm, matrix->nrows, update_s, update_p, verbose);
             if(verbose){
                     std::cout << *cpr_solve_ptr << std::endl;
             }
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


void solve_regular(int n, const mwIndex * cols, mwIndex const * rows, const double * entries, 
        const std::vector<double> & b, std::vector<double> & x, double tolerance,
		   int maxiter, int & iters, double & error, boost::property_tree::ptree prm){

    std::string relaxParam;
    /***************************************
     * Start AMGCL-link and select options *
     ***************************************/

    //boost::property_tree::ptree prm;
    /* Set tolerance */
    prm.put("solver.tol", tolerance);
    if(maxiter > 0){
        prm.put("solver.maxiter", maxiter);
    }
       
    /***************************************
     *        Solve problem                *
     ***************************************/
    auto t1 = std::chrono::high_resolution_clock::now();
    bool write_params = prm.get<bool>("write_params");
    if(write_params){
      std::cout << "Writing amgcl setup file to mrst_amgcl_cpr_setup.json xx" << std::endl;
      std::ofstream file("mrst_regular_setup.json");
      boost::property_tree::json_parser::write_json(file, prm);
    }
    int block_size = prm.get<int>("block_size");
    const auto matrix = amgcl::adapter::zero_copy(n, &cols[0], &rows[0], &entries[0]);
    bool verbose = prm.get<int>("verbosity") >0 ;
    switch(block_size){
      case 0:
      case 1:
      {
        auto t2 = std::chrono::high_resolution_clock::now();
        if(verbose){
            std::cout << "Solver setup took "
                      << (double)std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count()/1000.0
                      << " seconds\n";
        }
	std::tie(iters, error) = solve_shared(scalar_solve_ptr, matrix, b, x, prm, verbose);
        if(verbose){
            std::cout << *scalar_solve_ptr << std::endl;
        }
      } break;
       BOOST_PP_SEQ_FOR_EACH(AMGCL_BLOCK_SOLVER, block_solve_ptr, AMGCL_BLOCK_SIZES)
      default:
       {
	 std::cout << "Block size is :" << block_size << std::endl;
	 std::string msg("AMGCL:UndefBlockSize", "Failure: Block size not supported. xxx");	     
	 mexErrMsgIdAndTxt("AMGCL:UndefBlockSize",msg.c_str());
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
 

    double * entries;
    std::string relaxParam;
    std::string coarsenParam;

    if (nrhs != 5) {
      std::string msg("5 input arguments required. input is:");
      std::cout << "Given input " << nrhs << std::endl;
      mexErrMsgTxt(msg.c_str());
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
    double tolerance = mxGetScalar(prhs[2]);
    int maxiter = mxGetScalar(prhs[3]);
    boost::property_tree::ptree prm;
    {
      char *chopt = mxArrayToString(prhs[4]);  
      std::istringstream str(chopt);
      boost::property_tree::read_json(str, prm);
      std::cout << "options in prm" << std::endl;
      std::ofstream file("options_amgcl.json");
      //pt::write_json(std::cout, prm);
      boost::property_tree::write_json(file, prm);
    }

    std::vector<double> b(n);
    #pragma omp parallel for
    for(int ix = 0; ix < n; ix++){
        b[ix] = rhs[ix];
    }
    int    iters;
    double error;
    std::vector<double> x(M, 0.0);
    std::string solver_type = prm.get<std::string>("solver_type");
    int reuse_mode = prm.get<int>("reuse_mode");
    bool verbose = prm.get<int>("verbosity") >0 ;
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

    if( solver_type == "regular"){
      solve_regular(M, cols, rows, entries, b, x, tolerance, maxiter, iters, error, prm);
    }else if(solver_type == "cpr"){
      solve_cpr(M, cols, rows, entries, b, x, tolerance, maxiter, iters, error, prm);
    }else if(solver_type == "reset"){
	// Remove shared pointers
	if(verbose){
	    std::cout << "Resetting all solvers." << std::endl;
	}
	reset_solvers();
    }else{
      std::string msg("Unknown solver_type ");
      msg += solver_type;
      mexErrMsgTxt(msg.c_str());
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
    return;
}
