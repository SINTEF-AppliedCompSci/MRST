#include <fstream> //@@ for debugging
#include <iostream>
#include <string>
#include <algorithm>
#include <fstream>

#include "TensorComp.hpp"
#include "binaryop_algo.hpp"
#include "mex.hpp"
#include "mexAdapter.hpp"

using namespace matlab::data;
using matlab::mex::ArgumentList;

using namespace std;

class MexFunction : public matlab::mex::Function {

public:
  MexFunction() : matlabPtr(getEngine()) {}

  void operator() (ArgumentList outputs, ArgumentList inputs) {

    const SparseArray<double> sarr = inputs[0];

    for (auto it = sarr.begin(); it != sarr.end(); ++it) {
      SparseIndex index = sarr.getIndex(it);
      cout << "[" << index.first << ", " << index.second << "] - " << *it << '\n';
    }
  
    // create result array
    
    vector<double> elements(sarr.begin(), sarr.end());
    for_each(elements.begin(), elements.end(), [](double&d) {d *= 2;});
    size_t nnz = elements.size();

    vector<size_t> rows, cols;
    for (auto it = sarr.begin(); it != sarr.end(); ++it){
      rows.push_back(sarr.getIndex(it).first);
      cols.push_back(sarr.getIndex(it).second);
    }
      
    
    ArrayFactory factory;
    auto data_p = factory.createBuffer<double>(nnz);
    auto rows_p = factory.createBuffer<size_t>(nnz);
    auto cols_p = factory.createBuffer<size_t>(nnz);

    copy(elements.begin(), elements.end(), data_p.get());
    copy(rows.begin(), rows.end(), rows_p.get());
    copy(cols.begin(), cols.end(), cols_p.get());

    const ArrayDimensions adim = sarr.getDimensions();
    
    SparseArray<double> result =
      factory.createSparseArray<double> ( adim, nnz, move(data_p), move(rows_p), move(cols_p));

    outputs[0] = result;
    
  }
  
  // SparseArray result = factory.createSparseArray(

  
  // void operator() (ArgumentList outputs, ArgumentList inputs) {

  //   const Array arr = inputs[0];
  //   const ArrayDimensions adim = arr.getDimensions();

  //   auto elements = getReadOnlyElements<double>(arr);
  //   auto it = elements.begin();
    
  //   for (int rows = 0; rows != adim[0]; ++rows) {
  //     for (int cols = 0; cols != adim[1]; ++cols) {
  //       cout << *it++ << " ";
  //     }
  //     cout << endl;
  //   }
  // }
  
private:
  std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr;
};
