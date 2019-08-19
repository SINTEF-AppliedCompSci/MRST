template <class T, typename M>
std::tuple<size_t, double> solve_shared(std::shared_ptr<T> & solve_ptr,
                          const M matrix,
                          const std::vector<double> & b,
                          std::vector<double> & x,
                          boost::property_tree::ptree & prm,
                          bool verbose){
      auto t1 = std::chrono::high_resolution_clock::now();
      bool initialized = solve_ptr != nullptr;
      bool do_setup;

      std::cout << matrix->nrows << std::endl;
      if(initialized){
        do_setup = check_preconditioner_validity(solve_ptr, matrix);
      }else{
        do_setup = true;
      }

      if(do_setup){
        std::cout << "Initializing solver!" << std::endl;
        solve_ptr = std::make_shared<T>(*matrix, prm);
      }
      auto t2 = std::chrono::high_resolution_clock::now();
      if(verbose){
          std::cout << "Solver setup took "
                    << (double)std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count()/1000.0
                    << " seconds\n";
      }
      std::tuple<size_t, double> result  = (*solve_ptr)(b, x);
      // std::tie(iters, error) = (*solve_ptr)(b, x);

      if(verbose){
          std::cout << (*solve_ptr) << std::endl;
      }
      return result;
};
template <class T, typename M>
bool check_preconditioner_validity(std::shared_ptr<T> & solve_ptr, const M matrix){
  size_t nrows_precond = (*solve_ptr).system_matrix_ptr()->nrows;
  size_t nrows = matrix -> nrows;
  std::cout << "Matrix: " << nrows << " Precond: " << nrows_precond << std::endl;
  // std::cout << (*solve_ptr).system_matrix_ptr()->nrows;
  // do_setup = matrix->nrows ==
  bool do_setup = nrows_precond != nrows;
  return do_setup;
}
