#include "mrst_duneistl.hpp"
#include <iostream>
int main()
{
  std::cout << "Hello, World!";
  constexpr int bz=1;
  mrst::BlockIlu0Solver<bz> solver;
  std::string matrixfile("matrix.txt");
  std::string rhsfile("rhs.txt");
  double* result;
  solver.solve(result,matrixfile, rhsfile);
  return 0;
}
