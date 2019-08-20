#define AMGCL_BLOCK_SOLVER(z, data, B)                                          \
  case B:                                                                       \
  {                                                                             \
  std::tie(iters, error) =                                                      \
  solve_shared(BOOST_PP_CAT(data, B), matrix, b, x, prm, verbose);              \
} break;
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
  static std::shared_ptr<ScalarSolver> BOOST_PP_CAT(block_solve_ptr, B)(nullptr);

// Define reset of block solver
#define AMGCL_RESET_BLOCK_SOLVER(z, data, B) BOOST_PP_CAT(block_solve_ptr, B).reset();
