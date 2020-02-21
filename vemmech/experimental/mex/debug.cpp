#include <fstream> //@@ for debugging
#include <iostream>
#include <string>
#include <algorithm>

#include "TensorComp.hpp"
#include "contract_algo.hpp"

using namespace std;

int main()
{
  ifstream is("saved.comps");

  int num_comps;
  is >> num_comps;

  vector<TensorComp<double>> comps(num_comps);
  for (int i = 0; i != num_comps; ++i)
    comps[i].read(is);

  const vector<TensorComp<double>> resultcomps = contract_components(comps);
  
}
