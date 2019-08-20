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

      if(initialized){
        do_setup = check_preconditioner_validity(solve_ptr, matrix, matrix->nrows);
        if(verbose && do_setup){
          std::cout << "Attempted to reuse preconditioner, but dimensions did not match." << std::endl;
        }
      }else{
        do_setup = true;
      }

      if(do_setup){
        if(verbose){
          std::cout << "Initializing solver..." << std::endl;
        }
        solve_ptr = std::make_shared<T>(*matrix, prm);
      }
      auto t2 = std::chrono::high_resolution_clock::now();
      if(verbose){
          std::cout << "Solver setup took "
                    << (double)std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count()/1000.0
                    << " seconds\n";
      }
      std::tuple<size_t, double> result  = (*solve_ptr)(b, x);

      if(verbose){
          std::cout << (*solve_ptr) << std::endl;
      }
      return result;
};

template <class T, typename M, typename V, typename W>
std::tuple<size_t, double> solve_shared_cpr(std::shared_ptr<T> & solve_ptr,
                          const M matrix,
                          const V & b,
                          W & x,
                          boost::property_tree::ptree & prm,
                          size_t nrows,
                          bool update_sprecond,
                          bool verbose){
      auto t1 = std::chrono::high_resolution_clock::now();
      bool initialized = solve_ptr != nullptr;
      bool do_setup;
      if(initialized){
        do_setup = check_preconditioner_validity(solve_ptr, matrix, nrows);
        if(!do_setup && update_sprecond){
          if(verbose){
            std::cout << "Updating second-stage preconditioner only." << std::endl;
            auto cpr = solve_ptr->precond();
            cpr.update_sprecond(matrix);
          }
        }
        if(verbose && do_setup){
          std::cout << "Attempted to reuse preconditioner, but dimensions did not match." << std::endl;
        }
      }else{
        do_setup = true;
      }

      if(do_setup){
        if(verbose){
          std::cout << "Initializing CPR..." << std::endl;
        }
        solve_ptr = std::make_shared<T>(matrix, prm);
      }
      auto t2 = std::chrono::high_resolution_clock::now();
      if(verbose){
          std::cout << "Solver setup took "
                    << (double)std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count()/1000.0
                    << " seconds\n";
      }
      std::tuple<size_t, double> result  = (*solve_ptr)(b, x);

      if(verbose){
          std::cout << (*solve_ptr) << std::endl;
      }
      return result;
};

template <class T, typename M>
bool check_preconditioner_validity(std::shared_ptr<T> & solve_ptr, const M matrix, size_t nrows){
  size_t nrows_precond = solve_ptr->system_matrix_ptr()->nrows;
  bool do_setup = nrows_precond != nrows;
  return do_setup;
}
