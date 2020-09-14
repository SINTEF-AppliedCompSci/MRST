#include <fstream> //@@ for debugging
#include <iostream>
#include <string>
#include <algorithm>

#include "TensorComp.hpp"
#include "binaryop_algo.hpp"
#include "adi_helpers.hpp"

#include "mex.hpp"
#include "mexAdapter.hpp"

using namespace matlab::data;
using matlab::mex::ArgumentList;

using namespace std;

class MexFunction : public matlab::mex::Function {

public:
  MexFunction() : matlabPtr(getEngine()) {}

  // Input should be:
  // inputs[0] - cell array with two tensor components that are /compatible/,
  //             i.e. they have the same index names and index orders.
  // inputs[1] - string indicating the binary operation ('-', '+', '*', '/')
  void operator() (ArgumentList outputs, ArgumentList inputs) {

    checkArguments(outputs, inputs);

    const CellArray cells = inputs[0];
    const CharArray op = inputs[1];
    if (cells.getNumberOfElements() != 2)
      raiseError("tbinaryop takes exactly two components as input.");
    
    // if either of the two components uses ADI, then ADI must be used throughout
    const bool use_adi = check_adi(cells[0]) || check_adi(cells[1]);

    if (use_adi)
      outputs[0] = apply_operator<BasicAD>(cells, op);
    else
      outputs[0] = apply_operator<double>(cells, op);

    //std::cout << "Finished binary operator - now returning to MATLAB." << std::endl;
  }

  template<typename T>
  StructArray apply_operator(const CellArray& cells, const CharArray& op) {

    const auto c1 = TensorComp<T>(extract_indexnames(cells[0]),
                                  extract_numbers<T>(cells[0], "coefs", matlabPtr),
                                  extract_numbers<size_t>(cells[0], "ixs", matlabPtr));
    const auto c2 = TensorComp<T>(extract_indexnames(cells[1]),
                                  extract_numbers<T>(cells[1], "coefs", matlabPtr),
                                  extract_numbers<size_t>(cells[1], "ixs", matlabPtr));

    const TensorComp<T> resultcomp =
      apply_binary_op(c1, c2, string(op.begin(), op.end()));

    // now convert TensorComp to a return value for MATLAB
    //cout << "applied binary_op, now preparing data for returning to MATLAB." << endl;
    ArrayFactory factory;
    StructArray result = factory.createStructArray({1,1},
                                                   {"indexnames", "coefs", "ixs"});
    CellArray indexnames = factory.createCellArray({1, resultcomp.numIndices()});

    for (int j = 0; j != resultcomp.numIndices(); ++j)
      indexnames[j] = factory.createCharArray(resultcomp.indexNames()[j]);

    
    Array coefs = extract_coefs(resultcomp, matlabPtr);

    // TypedArray<double> coefs =
    //   factory.createArray<double>(ArrayDimensions {resultcomp.coefs().size(), 1},
    //                               &resultcomp.coefs()[0],
    //                               &resultcomp.coefs()[0] + resultcomp.coefs().size());

    
    vector<double> ixs_double(resultcomp.ixs().begin(), resultcomp.ixs().end());
    TypedArray<double> ixs =
      factory.createArray<double>(
        ArrayDimensions {resultcomp.coefs().size(), resultcomp.numIndices()},
        &ixs_double[0], &ixs_double[0] + resultcomp.ixs().size());

    result[0]["indexnames"] = indexnames;
    result[0]["coefs"] = coefs;
    result[0]["ixs"] = ixs;

    return result;
  }
  

    // ofstream os("saved.comps");
    // os << 2 << '\n';
    // c1.write(os);
    // c2.write(os);
    // os.close();

  // // ----------------------------------------------------------------------------
  // vector<string> extract_indexnames_(const StructArray& comp)  {
  //   const CellArray ixnames = comp[0][string("indexnames")];

  //   vector<string> result;    
  //   for (const CharArray name : ixnames) 
  //     result.emplace_back(string(name.begin(), name.end()));
    
  //   return result;
  // }
  
  // // ----------------------------------------------------------------------------
  // template<typename T>
  // vector<T> extract_numbers_(const StructArray& comp, const string field) {
  //   const TypedArray<double> arr = comp[0][field];
  //   vector<T> result;
  //   copy(arr.begin(), arr.end(), back_inserter(result));
  //   return result;
  // }
  
  
  void checkArguments(ArgumentList outputs, ArgumentList inputs) {
  
    if (inputs.empty() || inputs[0].getType() != ArrayType::CELL)
      raiseError("First argument should be a cell array with two components.");

    if (inputs[1].getType() != ArrayType::CHAR)
      raiseError("Second argument should be a string.");

  }

  void raiseError(const string& str) const
  {
    ArrayFactory factory;
    matlabPtr->feval(u"error",
                     0,
                     vector<Array>({ factory.createScalar(str)}));
  }

private:
  std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr;
};

  
