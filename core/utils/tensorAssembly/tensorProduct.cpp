#include "mex.hpp"
#include "mexAdapter.hpp"
#include <typeinfo>
#include <iostream>

using namespace matlab::data;
using matlab::mex::ArgumentList;

class MexFunction : public matlab::mex::Function {
public:
  void operator()(ArgumentList outputs, ArgumentList inputs) {

    // input[0] : integer, size of output
    // input[1] : vector of double, first input vector
    // input[2] : vector of double, second input vector
    // input[3] : vector of uint64, dispatch indices for first input vector
    // input[4] : vector of uint64, dispatch indices for second input vector
    // input[5] : vector of uint64, dispatch indices for output vector
    // Validate arguments
    
    checkArguments(outputs, inputs);

    uint n = inputs[0][0];

    const TypedArray<double>&  u1 = inputs[1];
    const TypedArray<double>&  u2 = inputs[2];

    const TypedArray<uint64_t>&  dispind1 = inputs[3];
    const TypedArray<uint64_t>&  dispind2 = inputs[4];
    const TypedArray<uint64_t>&  dispind3 = inputs[5];
    
    ArrayDimensions dims = {n};
    
    ArrayFactory factory;
    
    TypedArray<double> v = factory.createArray<double>(dims);

    
    for (size_t i = 0; i < dispind1.getNumberOfElements(); ++i) {
      
      v[dispind3[i] - 1] +=  u1[dispind1[i] - 1] * u2[dispind2[i] - 1];
      
    }

    // Assign outputs
    outputs[0] = v;
    
  }

  void checkArguments(ArgumentList outputs, ArgumentList inputs) {
      
    std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr = getEngine();
        
    ArrayFactory factory;

    for (size_t i = 1; i < 3; ++i) {
      if (inputs[i].getType() != ArrayType::DOUBLE) {
        matlabPtr->feval(u"error", 0, 
                         std::vector<Array>({ factory.createScalar("Input must be double array") }));
      }
    }

    for (size_t i = 3; i < 6; ++i) {
      if (inputs[i].getType() != ArrayType::UINT64) {
        matlabPtr->feval(u"error", 0, 
                         std::vector<Array>({ factory.createScalar("Input must be uint64 array") }));
      }
    }

    if (outputs.size() > 1) {
      matlabPtr->feval(u"error", 0, 
                       std::vector<Array>({ factory.createScalar("Only one output is returned") }));
    }
  }
};
