#include <fstream> //@@ for debugging
#include <assert.h> //@@ for debugging
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
  MexFunction() {}; //: matlabPtr(getEngine()) {}

  // Input should be:
  // inputs[0] - cell array with two tensor components that are /compatible/,
  //             i.e. they have the same index names and index orders.
  // inputs[1] - string indicating the binary operation ('-', '+', '*', '/')
  void operator() (ArgumentList outputs, ArgumentList inputs) {

    //std::cout << "Check arguments." << std::endl;
    checkArguments(outputs, inputs);
    
    const CellArray cells = inputs[0];
    const CharArray op = inputs[1];

    //std::cout << "getNumberOfELements." << std::endl;
    if (cells.getNumberOfElements() != 2)
      raiseError("tbinaryop takes exactly two components as input.");

    // if either of the two components uses ADI, then ADI must be used throughout
    //std::cout << "check adi." << std::endl;
    const bool use_adi = check_adi(cells[0]) || check_adi(cells[1]);

    if (use_adi) {
      //std::cout << "apply adi" << std::endl;
      outputs[0] = apply_operator<BasicAD>(cells, op);
    } else { 
      //std::cout << "apply double" << std::endl;
      outputs[0] = apply_operator<double>(cells, op);
    }

    //std::cout << "Finished binary operator - now returning to MATLAB." << std::endl;
  }

  template<typename T>
  StructArray apply_operator(const CellArray& cells, const CharArray& op) {

    //std::cout << "1" << std::endl;
    const auto c1 = TensorComp<T>(extract_indexnames(cells[0]),
                                  extract_numbers<T>(cells[0], "coefs", getEngine()),
                                  extract_numbers<size_t>(cells[0], "ixs", getEngine()));
    //std::cout << "2" << std::endl;
    const auto c2 = TensorComp<T>(extract_indexnames(cells[1]),
                                  extract_numbers<T>(cells[1], "coefs", getEngine()),
                                  extract_numbers<size_t>(cells[1], "ixs", getEngine()));
    //std::cout << "3" << std::endl;

    const TensorComp<T> resultcomp =
      apply_binary_op(c1, c2, string(op.begin(), op.end()));

    //std::cout << "4" << std::endl;
    // now convert TensorComp to a return value for MATLAB
    //cout << "applied binary_op, now preparing data for returning to MATLAB." << endl;
    ArrayFactory factory;
    StructArray result = factory.createStructArray({1,1},
                                                   {"indexnames", "coefs", "ixs"});
    CellArray indexnames = factory.createCellArray({1, resultcomp.numIndices()});

    for (int j = 0; j != resultcomp.numIndices(); ++j)
      indexnames[j] = factory.createCharArray(resultcomp.indexNames()[j]);

    
    Array coefs = extract_coefs(resultcomp, getEngine());

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
    //std::cout << "5" << std::endl;
    
    return result;
  }
  
  void checkArguments(ArgumentList outputs, ArgumentList inputs) {

    if (inputs.size() != 2)
      raiseError("There should be exactly two arguments to tbinaryop.");
    
    if (inputs.empty() || inputs[0].getType() != ArrayType::CELL)
      raiseError("First argument should be a cell array with two components.");

    if (inputs[1].getType() != ArrayType::CHAR)
      raiseError("Second argument should be a string.");
  }

  void raiseError(const string& str) 
  {
    ArrayFactory factory;
    getEngine()->feval(u"error",
                       0,
                       vector<Array>({ factory.createScalar(str)}));
  }

private:
  // std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr;
};

  
