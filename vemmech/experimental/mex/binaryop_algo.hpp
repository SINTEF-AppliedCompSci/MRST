#ifndef __BINARYOP_ALGO_HPP
#define __BINARYOP_ALGO_HPP

#include <algorithm>
#include <string>
#include <vector>
#include <stdexcept>
#include <assert.h>

#include "TensorComp.hpp"
#include "multiix_tools.hpp"



// ----------------------------------------------------------------------------
template<typename It> bool ixless(const It& iter1, const It& iter2, int N)
// ----------------------------------------------------------------------------
{
  for (int i = 0; i != N; ++i)
    if (iter1[i] != iter2[i])
      return (iter1[i] < iter2[i]);
  return false;
}

// ----------------------------------------------------------------------------
// this function requires t1 and t2 to be already sorted.  It should only be
// called through a dispatch from the 'apply_binary_op' function.
template<typename T> inline TensorComp<T>
apply_add_op(const TensorComp<T>& t1, const TensorComp<T>& t2)
// ----------------------------------------------------------------------------
{
  typedef typename TensorComp<T>::Index Index;
  const int Ni = t1.numIndices();
  assert(t2.numIndices() == Ni); // components should be compatible
  
  std::vector<Index> ixs1(t1.ixs());  transpose(ixs1, Ni);  // @@ perhaps reimplement if
  std::vector<Index> ixs2(t2.ixs());  transpose(ixs2, Ni);  // copying become too memory-intensive 

  assert(ixs1.size() == Ni * t1.coefs().size());
  assert(ixs2.size() == Ni * t2.coefs().size());
  
  auto t1_ix_iter = ixs1.begin();
  auto t1_ix_iter_end = ixs1.end();
  auto t1_coef_iter = t1.coefs().begin();

  auto t2_ix_iter = ixs2.begin();
  auto t2_ix_iter_end = ixs2.end();
  auto t2_coef_iter = t2.coefs().begin();

  std::vector<T> target_coefs;
  std::vector<typename TensorComp<T>::Index> target_ixs;
  
  while (t1_ix_iter != t1_ix_iter_end ||
         t2_ix_iter != t2_ix_iter_end) {
    
    if (t1_ix_iter == t1_ix_iter_end) {
      target_coefs.push_back(*t2_coef_iter++);
      for (int i = 0; i != Ni; ++i)
        target_ixs.push_back(*t2_ix_iter++);
    } else if (t2_ix_iter == t2_ix_iter_end) {
      target_coefs.push_back(*t1_coef_iter++);
      for (int i = 0; i != Ni; ++i)
        target_ixs.push_back(*t1_ix_iter++);
    } else if (ixless(t1_ix_iter, t2_ix_iter, Ni)) {
      // add 0 to the coefficient from t1
      target_coefs.push_back(*t1_coef_iter++);
      for (int i = 0; i != Ni; ++i)
        target_ixs.push_back(*t1_ix_iter++);
    } else if (ixless(t2_ix_iter, t1_ix_iter, Ni)) {
      // add 0 to the coefficient from t2
      target_coefs.push_back(*t2_coef_iter++);
      for (int i = 0; i != Ni; ++i)
        target_ixs.push_back(*t2_ix_iter++);
    } else {
      // the two indices are equal
      target_coefs.push_back(*t1_coef_iter + *t2_coef_iter);
      for (int i = 0; i != Ni; ++i)
        target_ixs.push_back(*t1_ix_iter++);
      
      ++t1_coef_iter;
      ++t2_coef_iter;
      t2_ix_iter += Ni; // t2_ix_iter was already advanced in the loop above
    }
  }
  
  assert(target_coefs.size() * Ni == target_ixs.size());
  transpose(target_ixs, target_coefs.size());
  return TensorComp<T>(t1.indexNames(), target_coefs, target_ixs);
}

// ----------------------------------------------------------------------------
// this function requires t1 and t2 to be already sorted.  It should only be
// called through a dispatch from the 'apply_binary_op' function.
template<typename T> inline TensorComp<T>
apply_subtract_op(const TensorComp<T>& t1, const TensorComp<T>& t2)
// ----------------------------------------------------------------------------
{
  typedef typename TensorComp<T>::Index Index;
  const int Ni = t1.numIndices();
  assert(t2.numIndices() == Ni); // components should be compatible

  std::vector<Index> ixs1(t1.ixs());  transpose(ixs1, Ni);  // @@ perhaps reimplement if
  std::vector<Index> ixs2(t2.ixs());  transpose(ixs2, Ni);  // copying become too memory-intensive 

  assert(ixs1.size() == Ni * t1.coefs().size());
  assert(ixs2.size() == Ni * t2.coefs().size());
  
  auto t1_ix_iter = ixs1.begin();
  auto t1_ix_iter_end = ixs1.end();
  auto t1_coef_iter = t1.coefs().begin();

  auto t2_ix_iter = ixs2.begin();
  auto t2_ix_iter_end = ixs2.end();
  auto t2_coef_iter = t2.coefs().begin();

  std::vector<T> target_coefs;
  std::vector<typename TensorComp<T>::Index> target_ixs;
                      
  while (t1_ix_iter != t1_ix_iter_end ||
         t2_ix_iter != t2_ix_iter_end) {

    if (t1_ix_iter == t1_ix_iter_end) {
      target_coefs.push_back(-1 * (*t2_coef_iter++));
      for (int i = 0; i != Ni; ++i)
        target_ixs.push_back(*t2_ix_iter++);
    } else if (t2_ix_iter == t2_ix_iter_end) {
      target_coefs.push_back(*t1_coef_iter++);
      for (int i = 0; i != Ni; ++i)
        target_ixs.push_back(*t1_ix_iter++);
    } else if (ixless(t1_ix_iter, t2_ix_iter, Ni)) {
      // subtract 0
      target_coefs.push_back(*t1_coef_iter++);
      for (int i = 0; i != Ni; ++i)
        target_ixs.push_back(*t1_ix_iter++);
    } else if (ixless(t2_ix_iter, t1_ix_iter, Ni)) {
      // subtract t2 from 0
      target_coefs.push_back(-1 * (*t2_coef_iter++));
      for (int i = 0; i != Ni; ++i)
        target_ixs.push_back(*t2_ix_iter++);
    } else {
      // subtract t2 from t1
      target_coefs.push_back(*t1_coef_iter - *t2_coef_iter);
      for (int i = 0; i != Ni; ++i)
        target_ixs.push_back(*t1_ix_iter++);

      ++t1_coef_iter;
      ++t2_coef_iter;
      t2_ix_iter += Ni; // t2_ix_iter was already advanced in the loop above
    }
  }
  assert(target_coefs.size() == Ni * target_ixs.size());
  transpose(target_ixs, target_coefs.size());
  return TensorComp<T>(t1.indexNames(), target_coefs, target_ixs);
}

// ----------------------------------------------------------------------------
// this function requires t1 and t2 to be already sorted.  It should only be
// called through a dispatch from the 'apply_binary_op' function.
template<typename T> inline TensorComp<T>
apply_mul_op(const TensorComp<T>& t1, const TensorComp<T>& t2)
// ----------------------------------------------------------------------------
{
  typedef typename TensorComp<T>::Index Index;
  const int Ni = t1.numIndices();
  assert(t2.numIndices() == Ni); // components should be compatible

  std::vector<Index> ixs1(t1.ixs());  transpose(ixs1, Ni);  // @@ perhaps reimplement if
  std::vector<Index> ixs2(t2.ixs());  transpose(ixs2, Ni);  // copying become too memory-intensive 

  assert(ixs1.size() == Ni * t1.coefs().size());
  assert(ixs2.size() == Ni * t2.coefs().size());
  
  auto t1_ix_iter = ixs1.begin();
  auto t1_ix_iter_end = ixs1.end();
  auto t1_coef_iter = t1.coefs().begin();

  auto t2_ix_iter = ixs2.begin();
  auto t2_ix_iter_end = ixs2.end();
  auto t2_coef_iter = t2.coefs().begin();

  std::vector<T> target_coefs;
  std::vector<typename TensorComp<T>::Index> target_ixs;

  while (t1_ix_iter != t1_ix_iter_end &&
         t2_ix_iter != t2_ix_iter_end) {

    if (ixless(t1_ix_iter, t2_ix_iter, Ni)) {
      // multiplication by 0 -> no new coefficient, but advance index
      t1_ix_iter += Ni;
      ++t1_coef_iter;
    } else if (ixless(t2_ix_iter, t1_ix_iter, Ni)) {
      // multiplication by 0 -> no new coefficient, but advance index
      t2_ix_iter += Ni;
      ++t2_coef_iter;
    } else {
      // the two indices are equal
      target_coefs.push_back(*t1_coef_iter * *t2_coef_iter);
      for (int i = 0; i != Ni; ++i)
        target_ixs.push_back(*t1_ix_iter++);

      ++t1_coef_iter;
      ++t2_coef_iter;
      t2_ix_iter += Ni; // t2_ix_iter was already advanced in the loop above
    }
  }

  assert(target_coefs.size() == Ni * target_ixs.size());
  
  transpose(target_ixs, target_coefs.size());
  return TensorComp<T>(t1.indexNames(), target_coefs, target_ixs);
}

// ----------------------------------------------------------------------------
// this function requires t1 and t2 to be already sorted.  It should only be
// called through a dispatch from the 'apply_binary_op' function.
template<typename T> inline TensorComp<T>
apply_div_op(const TensorComp<T>& t1, const TensorComp<T>& t2)
// ----------------------------------------------------------------------------
{
  typedef typename TensorComp<T>::Index Index;
  const int Ni = t1.numIndices();
  assert(t2.numIndices() == Ni); // components should be compatible

  std::vector<Index> ixs1(t1.ixs());  transpose(ixs1, Ni);  // @@ perhaps reimplement if
  std::vector<Index> ixs2(t2.ixs());  transpose(ixs2, Ni);  // copying become too memory-intensive 

  assert(ixs1.size() == Ni * t1.coefs().size());
  assert(ixs2.size() == Ni * t2.coefs().size());
  
  auto t1_ix_iter = ixs1.begin();
  auto t1_ix_iter_end = ixs1.end();
  auto t1_coef_iter = t1.coefs().begin();

  auto t2_ix_iter = ixs2.begin();
  auto t2_ix_iter_end = ixs2.end();
  auto t2_coef_iter = t2.coefs().begin();

  std::vector<T> target_coefs;
  std::vector<typename TensorComp<T>::Index> target_ixs;

  while (t1_ix_iter != t1_ix_iter_end &&
         t2_ix_iter != t2_ix_iter_end) {

    if (ixless(t1_ix_iter, t2_ix_iter, Ni)) {
      // division by zero
      throw std::runtime_error("Division by zero in element-wise tensor division.");
      
    } else if (ixless(t2_ix_iter, t1_ix_iter, Ni)) {
      // 0 divided by something -> no new coefficient, but advance index
      t2_ix_iter += Ni;
      ++t2_coef_iter;
    } else {
      // the two indices are equal
      target_coefs.push_back(*t1_coef_iter / *t2_coef_iter);
      for (int i = 0; i != Ni; ++i)
        target_ixs.push_back(*t1_ix_iter++);

      ++t1_coef_iter;
      ++t2_coef_iter;
      t2_ix_iter += Ni; // t2_ix_iter was already advanced in the loop above
    }
  }

  assert(target_coefs.size() == Ni * target_ixs.size());
  
  transpose(target_ixs, target_coefs.size());
  return TensorComp<T>(t1.indexNames(), target_coefs, target_ixs);
}

// ============================================================================
template<typename T> inline TensorComp<T>
apply_binary_op(TensorComp<T> t1,
                TensorComp<T> t2,
                const std::string operator_name)
// ============================================================================
{
  typedef typename TensorComp<T>::Index Index;

  t1.sortElementsByIndex();
  t2.sortElementsByIndex();

  assert(t1.numIndices() == t2.numIndices());
  
  // do operator dispatch at once, to avoid checks at each iteration
  // (e.g. whether the operator is symmetrical)
  if (operator_name == "+")
    return apply_add_op(t1, t2);
  else if (operator_name == "-")
    return apply_subtract_op(t1, t2);
  else if (operator_name == "*")
    return apply_mul_op(t1, t2);
  else if (operator_name == "/")
    return apply_div_op(t1, t2);
  else
    throw std::runtime_error("Wrong operator name, must be '+', '-', '*' or '/'.");

  return t1;
}

#endif
