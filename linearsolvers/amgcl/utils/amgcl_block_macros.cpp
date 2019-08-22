
// Define block matrix, block vector and corresponding block backends
#define AMGCL_DEFINE_BLOCK_TYPES(z, data, B)                                                        \
  typedef amgcl::static_matrix<double, B, B> BOOST_PP_CAT(BlockMat, B);                             \
  typedef amgcl::static_matrix<double, B, 1> BOOST_PP_CAT(BlockVec, B);                             \
  typedef amgcl::backend::builtin<BOOST_PP_CAT(BlockMat, B)> BOOST_PP_CAT(BlockBackend, B);
// Define block solver with preconditioner
#define AMGCL_DEFINE_BLOCK_SOLVER(z, data, B)                                                       \
  typedef amgcl::make_block_solver<                                                                 \
      amgcl::runtime::preconditioner<BOOST_PP_CAT(BlockBackend, B)>,                                \
      amgcl::runtime::solver::wrapper<BOOST_PP_CAT(BlockBackend, B)>                                \
  > BOOST_PP_CAT(data, B);                                                                          \
  static std::shared_ptr<BOOST_PP_CAT(data, B)> BOOST_PP_CAT(block_solve_ptr, B)(nullptr);

// Insert block solvers in switch
#define AMGCL_BLOCK_SOLVER(z, data, B)                                          \
  case B:                                                                       \
  {                                                                             \
  std::tie(iters, error) =                                                      \
  solve_shared(BOOST_PP_CAT(data, B), matrix, b, x, prm, verbose);              \
} break;

// Define block CPR
#define AMGCL_DEFINE_BLOCK_CPR_SOLVERS(z, data, B)                             \
  typedef amgcl::relaxation::as_preconditioner<                                \
        BOOST_PP_CAT(BlockBackend, B),                                         \
        amgcl::runtime::relaxation::wrapper                                    \
        >                                                                      \
      BOOST_PP_CAT(SPrecond, B);                                               \
  typedef amgcl::make_solver<                                                  \
      amgcl::preconditioner::cpr<PPrecond, BOOST_PP_CAT(SPrecond, B)>,         \
      amgcl::runtime::solver::wrapper<BOOST_PP_CAT(BlockBackend, B)>           \
      > BOOST_PP_CAT(CPRSolverBlock, B);                                       \
      static std::shared_ptr<BOOST_PP_CAT(CPRSolverBlock, B)>                  \
                        BOOST_PP_CAT(cpr_block_solve_ptr, B)(nullptr);         \
  typedef amgcl::make_solver<                                                  \
      amgcl::preconditioner::cpr_drs<PPrecond, BOOST_PP_CAT(SPrecond, B)>,     \
      amgcl::runtime::solver::wrapper<BOOST_PP_CAT(BlockBackend, B)>           \
      > BOOST_PP_CAT(CPR_DRSSolverBlock, B);                                   \
      static std::shared_ptr<BOOST_PP_CAT(CPR_DRSSolverBlock, B)>              \
                  BOOST_PP_CAT(cpr_drs_block_solve_ptr, B)(nullptr);

// Insert block solvers in switch
#define AMGCL_BLOCK_CPR_SOLVER(z, solver_name, B)                               \
  case B:                                                                       \
  {                                                                             \
  typedef BOOST_PP_CAT(BlockVec, B) bvec;                                       \
  typedef BOOST_PP_CAT(BlockMat, B) bmat;                                       \
  size_t n = matrix->nrows / B;                                                 \
  auto BM = amgcl::adapter::block_matrix<bmat>(*matrix);                         \
  std::vector<bvec> x_local(b.size(), amgcl::math::zero<bvec>());               \
  auto b_ptr = reinterpret_cast<const bvec*>(b.data());                         \
  auto b_local = amgcl::make_iterator_range(b_ptr, b_ptr + n);                  \
  std::tie(iters, error) =                                                      \
  solve_shared_cpr(BOOST_PP_CAT(solver_name, B), BM, b_local, x_local, prm, n, update_s, update_p, verbose); \
  auto x_data = x_local.data();                                                 \
  for(int i = 0; i < matrix->nrows; i++){                                       \
    x[i] = x_data[i / B](i % B);                                                \
  }                                                                             \
} break;

// Define reset of named block solver
  #define AMGCL_RESET_BLOCK_SOLVER(z, data, B) BOOST_PP_CAT(data, B).reset();
