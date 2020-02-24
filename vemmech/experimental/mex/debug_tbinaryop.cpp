#include <fstream> //@@ for debugging
#include <iostream>
#include <string>
#include <algorithm>

#include "TensorComp.hpp"
#include "binaryop_algo.hpp"

using namespace std;

int main()
{
  ifstream is("saved.comps");
  int num_comps;
  is >> num_comps; // should be 2

  TensorComp<double> c1, c2;
  c1.read(is);
  c2.read(is);

  const auto result = apply_binary_op(c1, c2, "+");

}
