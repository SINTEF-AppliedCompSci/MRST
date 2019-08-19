#define AMGCL_BLOCK_SOLVER(z, data, B)                                          \
  case B:                                                                       \
  {                                                                             \
  std::tie(iters, error) =                                                      \
  solve_shared(BOOST_PP_CAT(data, B), matrix, b, x, prm, verbose);              \
} break;

#define AMGCL_DEFINE_BLOCK_SOLVER(z, data, B)                                                       \
  typedef amgcl::make_block_solver<                                                                 \
      amgcl::runtime::preconditioner <amgcl::backend::builtin<amgcl::static_matrix<double, B, B>>>, \
      amgcl::runtime::solver::wrapper<amgcl::backend::builtin<amgcl::static_matrix<double, B, B>>>  \
  > BOOST_PP_CAT(data, B);                                                                          \
  static std::shared_ptr<ScalarSolver> BOOST_PP_CAT(block_solve_ptr, B)(nullptr);

#define AMGCL_RESET_BLOCK_SOLVER(z, data, B) BOOST_PP_CAT(block_solve_ptr, B).reset();
