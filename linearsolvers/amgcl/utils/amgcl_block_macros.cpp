#define AMGCL_BLOCK_SOLVER(z, data, B)                                           \
  case B:                                                                        \
  {                                                                              \
  typedef amgcl::backend::builtin<amgcl::static_matrix<double, B, B> > BBackend; \
  amgcl::make_block_solver<                                                      \
      amgcl::runtime::preconditioner<BBackend>,                                  \
      amgcl::runtime::solver::wrapper<BBackend>                                  \
  > solve(*matrix, prm);                                                   \
  auto t2 = std::chrono::high_resolution_clock::now();                           \
  if(verbose){                                                                   \
      std::cout << "Solver setup took "                                          \
                << (double)std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count()/1000.0 \
                << " seconds\n";                                                 \
  }                                                                              \
  std::tie(iters, error) = solve(b, x);                                          \
  if(verbose){                                                                   \
      std::cout << solve << std::endl;                                           \
  }                                                                              \
} break;

#define AMGCL_DEFINE_BLOCK_SOLVER(z, data, B)                                                       \
  typedef amgcl::make_block_solver<                                                                 \
      amgcl::runtime::preconditioner <amgcl::backend::builtin<amgcl::static_matrix<double, B, B>>>, \
      amgcl::runtime::solver::wrapper<amgcl::backend::builtin<amgcl::static_matrix<double, B, B>>>  \
  > BOOST_PP_CAT(data, B);                                                                          \
  static std::shared_ptr<ScalarSolver> BOOST_PP_CAT(block_solve_ptr, B)(nullptr);

#define AMGCL_RESET_BLOCK_SOLVER(z, data, B) BOOST_PP_CAT(block_solve_ptr, B).reset();
