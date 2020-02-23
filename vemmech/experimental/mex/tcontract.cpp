#include <fstream> //@@ for debugging
#include <iostream>
#include <string>
#include <algorithm>

#include "TensorComp.hpp"
#include "contract_algo.hpp"
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
    vector<TensorComp<double>> comps;

    for (const auto c : cells)
      comps.emplace_back( TensorComp<double>{extract_indexnames_(c),
                                             extract_numbers_<double>(c, "coefs"),
                                             extract_numbers_<int>(c, "ixs")});

    cout << "Number of components read: " << comps.size() << endl;

    // // saving components
    // ofstream os("saved.comps");
    // os << comps.size() << '\n';
    // for (const auto c : comps)
    //   c.write(os);
    // os.close();
    
    
    const vector<TensorComp<double>> resultcomps = contract_components(comps);

    // convert TensorComp to a return value to put into 'outputs'
    ArrayFactory factory;

    CellArray result = factory.createCellArray({1, resultcomps.size()});

    for (int i = 0; i != resultcomps.size(); ++i) {
      const auto& res = resultcomps[i];
      
      StructArray entry = factory.createStructArray({1,1},
                                                    {"indexnames", "coefs", "ixs"});
      
      CellArray indexnames = factory.createCellArray({1, res.numIndices()});
      for (int j = 0; j != res.numIndices(); ++j)
        indexnames[j] = factory.createCharArray(res.indexNames()[j]);

      TypedArray<double> coefs =
        factory.createArray<double>(ArrayDimensions{res.coefs().size(), 1},
                                    &res.coefs()[0],
                                    &res.coefs()[0] + res.coefs().size());

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

    outputs[0] = result;

  }

  // ----------------------------------------------------------------------------
  vector<string> extract_indexnames_(const StructArray& comp)  {
    const CellArray ixnames = comp[0][string("indexnames")];

    vector<string> result;    
    for (const CharArray name : ixnames) 
      result.emplace_back(string(name.begin(), name.end()));
    
    return result;
  }
  
  // ----------------------------------------------------------------------------
  template<typename T>
  vector<T> extract_numbers_(const StructArray& comp, const string field) {
    const TypedArray<double> arr = comp[0][field];

    vector<T> result;
    copy(arr.begin(), arr.end(), back_inserter(result));
    return result;
  }
  
  
  void checkArguments(ArgumentList outputs, ArgumentList inputs) {
  
    if (inputs.empty() || inputs[0].getType() != ArrayType::CELL)
      raiseError("This function takes a cell array.");

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

  
