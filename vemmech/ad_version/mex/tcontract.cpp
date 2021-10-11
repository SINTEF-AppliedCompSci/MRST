#include <fstream> //@@ for debugging
#include <iostream>
#include <string>
#include <algorithm>
#include <stdexcept>

#include "TensorComp.hpp"
#include "contract_algo.hpp"
#include "BasicAD.hpp"
#include "adi_helpers.hpp"

#include "mex.hpp"
#include "mexAdapter.hpp"


using namespace matlab::data;
using matlab::mex::ArgumentList;

using namespace std;

class MexFunction : public matlab::mex::Function {

public:
  MexFunction() : matlabPtr(getEngine()) {}
  
  void operator() (ArgumentList outputs, ArgumentList inputs) {

    checkArguments(outputs, inputs);

    const CellArray cells = inputs[0];

    // if at least one compoent uses ADI, then ADI must be used throughout
    bool use_adi = false;
    for (const auto c : cells)
      use_adi = use_adi || check_adi(c);

    if (use_adi) 
      outputs[0] = contract<BasicAD>(cells);
    else 
      outputs[0] = contract<double>(cells);
  }

  template<typename T>
  CellArray contract(const CellArray& cells)
  {
    vector<TensorComp<T>> comps;
    for (const auto c: cells) 
      comps.emplace_back( TensorComp<T>{ extract_indexnames(c),
                                          extract_numbers<T>(c, "coefs", matlabPtr),
                                          extract_numbers<size_t>(c, "ixs", matlabPtr)});
    const vector<TensorComp<T>> resultcomps = contract_components(comps);

    // convert TensorComp to a return value to put into 'outputs'
    ArrayFactory factory;
    CellArray result = factory.createCellArray({1, resultcomps.size()});

    for (int i = 0; i != resultcomps.size(); ++i) {
    
      const auto& res = resultcomps[i];
      StructArray entry = factory.createStructArray({1,1},
                                                    {"indexnames", "coefs", "ixs"});
      CellArray indexnames = factory.createCellArray({1, res.numIndices()});
      for (int j = 0; j != res.numIndices(); ++j) {
        indexnames[j] = factory.createCharArray(res.indexNames()[j]);
      }
      Array coefs = extract_coefs(res, matlabPtr);
      vector<double> ixs_double(res.ixs().begin(), res.ixs().end());
      TypedArray<double> ixs =
        factory.createArray<double>(
           ArrayDimensions{res.coefs().size(), res.numIndices()},
           &ixs_double[0], &ixs_double[0] + ixs_double.size());
      entry[0]["indexnames"] = indexnames;
      entry[0]["coefs"] = coefs;
      entry[0]["ixs"] = ixs;
      result[i] = entry;
    }
    return result;
  }


  
  // ----------------------------------------------------------------------------  
  void checkArguments(ArgumentList outputs, ArgumentList inputs) {
  
    if (inputs.empty() || inputs[0].getType() != ArrayType::CELL)
      raiseError("This function takes a cell array.");

  }
  
  // ----------------------------------------------------------------------------
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

