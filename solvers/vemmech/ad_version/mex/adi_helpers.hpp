#ifndef __ADI_HELPERS_HPP
#define __ADI_HELPERS_HPP

#include "mex.hpp"
#include "mexAdapter.hpp"

#include "BasicAD.hpp"

#include <fstream> // @@ debug
#include <assert.h> 
using namespace matlab::data;
using matlab::mex::ArgumentList;


// ----------------------------------------------------------------------------
bool check_adi(const StructArray& comp) {
  // @@ really only tests if it is an object handle
  assert(comp.getNumberOfElements() == 1);
  return
    comp[0]["coefs"].getType() == ArrayType::HANDLE_OBJECT_REF ||
    comp[0]["coefs"].getType() == ArrayType::VALUE_OBJECT;
      
}
  
// ----------------------------------------------------------------------------
inline std::vector<std::string> extract_indexnames(const StructArray& comp)  {
  assert(comp.getNumberOfElements() == 1);
  
  const CellArray ixnames = comp[0][std::string("indexnames")];
    
  std::vector<std::string> result;    
  for (const CharArray name : ixnames) 
    result.emplace_back(std::string(name.begin(), name.end()));
    
  return result;
}
  
// ----------------------------------------------------------------------------  
inline Array extract_coefs(const TensorComp<double>& res,
                           std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr) {
  assert(res.coefs().size() > 0); // @@ in general might be unnecessary, added for debug purposes
  ArrayFactory factory;
  TypedArray<double> coefs =
    factory.createArray<double>(ArrayDimensions{res.coefs().size(), 1},
                                &res.coefs()[0],
                                &res.coefs()[0] + res.coefs().size());
  return coefs;
}

// ----------------------------------------------------------------------------  
SparseArray<double> make_sparse_(const ArrayDimensions& adim,
                                 const std::vector<double>& elems,
                                 const std::vector<size_t>& rows,
                                 const std::vector<size_t>& cols)
{
  ArrayFactory factory;

  size_t nnz = elems.size();
  
  assert(nnz > 0); // @@ for debug
  assert(rows.size() == nnz);
  assert(cols.size() == nnz);
  
  auto data_p = factory.createBuffer<double>(nnz);
  auto rows_p = factory.createBuffer<size_t>(nnz);
  auto cols_p = factory.createBuffer<size_t>(nnz);

  copy(elems.begin(), elems.end(), data_p.get());
  copy(rows.begin(), rows.end(), rows_p.get());
  copy(cols.begin(), cols.end(), cols_p.get());

  SparseArray<double> result =
    factory.createSparseArray<double>(adim, nnz,
                                      move(data_p), move(rows_p), move(cols_p));
  return result;
}
  
// ----------------------------------------------------------------------------  
inline Array extract_coefs(const TensorComp<BasicAD>& res,
                           std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr) {

  ArrayFactory factory;
  const std::vector<BasicAD>& coefs = res.coefs();

  // exctracting values
  std::vector<double> vals;
  for (size_t i = 0; i != coefs.size(); ++i)
    vals.push_back(coefs[i].val());

  TypedArray<double> newcoefs =
    factory.createArray<double>(ArrayDimensions{coefs.size(), 1},
                                &vals[0], &vals[0] + vals.size());
  // extracting derivatives
  std::vector<double> elems;
  std::vector<size_t> rows, cols;
  size_t num_cols = 0;
  for (size_t i = 0; i != coefs.size(); ++i) {
    const auto& derivs = coefs[i].derivs();
    num_cols = (derivs.size() > num_cols) ? derivs.size() : num_cols;
    for (size_t j = 0; j != derivs.size(); ++j) {
      if (derivs[j] != 0) {
        elems.push_back(derivs[j]);
        rows.push_back(i);
        cols.push_back(j);
      }
    }
  }
  std::vector<Array> args;
  args.push_back(newcoefs);
  if (elems.size() == 0)
    args.push_back(factory.createScalar(num_cols)); // workaround memory issue with empty SparseArrays
  else
    args.push_back(make_sparse_(ArrayDimensions {coefs.size(), num_cols},
                                elems, rows, cols));
    
  Array result = matlabPtr->feval(u"makeadi", args);
  return result;
}

// ----------------------------------------------------------------------------  
inline void parse_sparse_(const Array& arr,
                          std::vector<double>& elems,
                          std::vector<size_t>& rows,
                          std::vector<size_t>& cols) {

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

// ----------------------------------------------------------------------------
template<typename T> inline
std::vector<T> extract_numbers(const StructArray& comp, const std::string field,
                          std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr) {
  assert(comp.getNumberOfElements() == 1);
  const TypedArray<double> arr = comp[0][field];
  std::vector<T> result;
  copy(arr.begin(), arr.end(), back_inserter(result));
  return result;
}

// ----------------------------------------------------------------------------  
template<> inline 
std::vector<BasicAD> extract_numbers(const StructArray& comp, const std::string field,
                                std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr) {
  // @@ unimplemented
  std::vector<BasicAD> result;
  assert(comp.getNumberOfElements() == 1);
  Array adinput = comp[0][field];

  // the input may either be double or ADI, but we will read it as ADI
  if (adinput.getType() == ArrayType::VALUE_OBJECT ||
      adinput.getType() == ArrayType::HANDLE_OBJECT_REF) { // input is ADI
      
    // read values and construct BasicAD (for the moment without derivatives)
    TypedArray<double> vals = matlabPtr->feval(u"getval", adinput);
    std::vector<double> vals_double(vals.begin(), vals.end());

    // read jacobian
    CellArray jac = matlabPtr->feval(u"getjac", adinput);
    if (jac.getNumberOfElements() != 1) 
      throw std::runtime_error("Currently only supports one jacobian matrix.");

    const ArrayDimensions dims = jac[0].getDimensions();
    const int num_derivs = dims[1];
    
    // create vector of BasicAD
    result.reserve(vals_double.size());
    for (size_t i = 0; i != vals_double.size(); ++i)
      result.push_back(BasicAD {vals_double[i], num_derivs});

    std::vector<double> jac_elems;
    std::vector<size_t> rows, cols;
    parse_sparse_(jac[0], jac_elems, rows, cols);

    for (size_t i = 0; i != jac_elems.size(); ++i) {
      assert(rows[i] < result.size());
      result[rows[i]].setDeriv(cols[i], jac_elems[i]);
    }
      
  } else {
    // input presumably double
    const TypedArray<double>& darr = adinput;
    for (auto it = darr.begin(); it != darr.end(); ++it)
      result.push_back(BasicAD {*it});
  } 
    
  return result;
}
  
#endif
