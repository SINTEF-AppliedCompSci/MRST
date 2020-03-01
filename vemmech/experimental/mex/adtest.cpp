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
    
    Array adinput(inputs[0]);
    cout << (int)adinput.getType() << endl;
    
    TypedArray<double> vals = matlabPtr->feval(u"getval", adinput);    
    //const TypedArray<double> vals = matlabPtr->getProperty(adinput, 6, u"val");
    
    vector<double> dvals(vals.begin(), vals.end());
    
    for (int i = 0; i != dvals.size(); ++i)
      cout << dvals[i] << " ";
    cout << endl;
    
    CellArray jac = matlabPtr->feval(u"getjac", adinput);
    
    cout << endl;
    cout << jac.getNumberOfElements() << endl;
    
    vector<double> elems;
    vector<size_t> rows, cols;
    parse_sparse_(jac[0], elems, rows, cols);
    
    for_each(elems.begin(), elems.end(), [](double& d) { d*=2; });
    
    // make sparse return array for return
    const ArrayDimensions adim = jac[0].getDimensions();
    SparseArray<double> sparr = make_sparse_(adim, elems, rows, cols);

    
    for_each(dvals.begin(), dvals.end(), [](double& d) { d*=2;});
    ArrayFactory factory;
    TypedArray<double> newcoefs =
      factory.createArray<double>(ArrayDimensions{dvals.size(), 1},
                                  &dvals[0], &dvals[0] + dvals.size());

    vector<Array> args;
    args.push_back(newcoefs);
    args.push_back(sparr);
    Array res = matlabPtr->feval(u"makeadi", args);
    //matlabPtr->feval(u"makeadi", newcoefs, sparr);
    outputs[0] = res;
    //  outputs[1] = sparr;
    
    
  }

  SparseArray<double> make_sparse_(const ArrayDimensions& adim,
                                   const vector<double>& elems,
                                   const vector<size_t>& rows,
                                   const vector<size_t>& cols)
  {
    ArrayFactory factory;
    auto data_p = factory.createBuffer<double>(elems.size());
    auto rows_p = factory.createBuffer<size_t>(elems.size());
    auto cols_p = factory.createBuffer<size_t>(elems.size());

    copy(elems.begin(), elems.end(), data_p.get());
    copy(rows.begin(), rows.end(), rows_p.get());
    copy(cols.begin(), cols.end(), cols_p.get());

    SparseArray<double> result =
      factory.createSparseArray<double>(adim, elems.size(),
                                        move(data_p), move(rows_p), move(cols_p));
    return result;
  }

  void parse_sparse_(const Array& arr,
                     vector<double>& elems,
                     vector<size_t>& rows,
                     vector<size_t>& cols) {

    const SparseArray<double> sarr = arr;

    elems.clear();
    rows.clear();
    cols.clear();
    
    for (auto it = sarr.begin(); it != sarr.end(); ++it) {
      SparseIndex index = sarr.getIndex(it);
      elems.push_back(*it);
      rows.push_back(index.first);
      cols.push_back(index.second);
      //cout << "[" << index.first << ", " << index.second << "] - " << *it << '\n';
    }
  }
  //   // ---------------------
    
  //   const SparseArray<double> sarr = inputs[0];

  //   for (auto it = sarr.begin(); it != sarr.end(); ++it) {
  //     SparseIndex index = sarr.getIndex(it);
  //     cout << "[" << index.first << ", " << index.second << "] - " << *it << '\n';
  //   }
  
  //   // create result array
    
  //   vector<double> elements(sarr.begin(), sarr.end());
  //   for_each(elements.begin(), elements.end(), [](double&d) {d *= 2;});
  //   size_t nnz = elements.size();

  //   vector<size_t> rows, cols;
  //   for (auto it = sarr.begin(); it != sarr.end(); ++it){
  //     rows.push_back(sarr.getIndex(it).first);
  //     cols.push_back(sarr.getIndex(it).second);
  //   }
      
    
  //   ArrayFactory factory;
  //   auto data_p = factory.createBuffer<double>(nnz);
  //   auto rows_p = factory.createBuffer<size_t>(nnz);
  //   auto cols_p = factory.createBuffer<size_t>(nnz);

  //   copy(elements.begin(), elements.end(), data_p.get());
  //   copy(rows.begin(), rows.end(), rows_p.get());
  //   copy(cols.begin(), cols.end(), cols_p.get());

  //   const ArrayDimensions adim = sarr.getDimensions();
    
  //   SparseArray<double> result =
  //     factory.createSparseArray<double> ( adim, nnz, move(data_p), move(rows_p), move(cols_p));

  //   outputs[0] = result;
    
  // }
  
  
private:
  std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr;
};
