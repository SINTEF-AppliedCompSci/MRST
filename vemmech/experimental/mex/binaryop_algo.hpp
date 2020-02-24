#ifndef __BINARYOP_ALGO_HPP
#define __BINARYOP_ALGO_HPP

#include <algorithm>
#include <string>
#include <vector>
#include <stdexcept>

#include "TensorComp.hpp"
#include "multiix_tools.hpp"

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
  
  auto t1_ix_iter = ixs1.begin();
  auto t1_ix_iter_end = ixs1.end();
  auto t1_coef_iter = t1.coefs().begin();

  auto t2_ix_iter = ixs2.begin();
  auto t2_ix_iter_end = ixs2.end();
  auto t2_coef_iter = t2.coefs().begin();

  std::vector<T> target_coefs;
  std::vector<typename TensorComp<T>::Index> target_ixs;

  const auto ixless =
    [Ni](const decltype(t1_ix_iter) iter1, const decltype(t2_ix_iter) iter2) {
      for (int i = 0; i != Ni; ++i)
        if (iter1[i] != iter2[i])
          return (iter1[i] < iter2[i]);
      return false;
    };
                      
  while (t1_ix_iter != t1_ix_iter_end &&
         t2_ix_iter != t2_ix_iter_end) {

    if (ixless(t1_ix_iter, t2_ix_iter)) {
      // add 0 to the coefficient from t1
      target_coefs.push_back(*t1_coef_iter++);
      for (int i = 0; i != Ni; ++i)
        target_ixs.push_back(*t1_ix_iter++);
    } else if (ixless(t2_ix_iter, t1_ix_iter)) {
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
  
  auto t1_ix_iter = ixs1.begin();
  auto t1_ix_iter_end = ixs1.end();
  auto t1_coef_iter = t1.coefs().begin();

  auto t2_ix_iter = ixs2.begin();
  auto t2_ix_iter_end = ixs2.end();
  auto t2_coef_iter = t2.coefs().begin();

  std::vector<T> target_coefs;
  std::vector<typename TensorComp<T>::Index> target_ixs;

  const auto ixless =
    [Ni](const decltype(t1_ix_iter) iter1, const decltype(t2_ix_iter) iter2) {
      for (int i = 0; i != Ni; ++i)
        if (iter1[i] != iter2[i])
          return (iter1[i] < iter2[i]);
      return false;
    };
                      
  while (t1_ix_iter != t1_ix_iter_end &&
         t2_ix_iter != t2_ix_iter_end) {

    if (ixless(t1_ix_iter, t2_ix_iter)) {
      // subtract 0
      target_coefs.push_back(*t1_coef_iter++);
      for (int i = 0; i != Ni; ++i)
        target_ixs.push_back(*t1_ix_iter++);
    } else if (ixless(t2_ix_iter, t1_ix_iter)) {
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
  
  auto t1_ix_iter = ixs1.begin();
  auto t1_ix_iter_end = ixs1.end();
  auto t1_coef_iter = t1.coefs().begin();

  auto t2_ix_iter = ixs2.begin();
  auto t2_ix_iter_end = ixs2.end();
  auto t2_coef_iter = t2.coefs().begin();

  std::vector<T> target_coefs;
  std::vector<typename TensorComp<T>::Index> target_ixs;

  const auto ixless =
    [Ni](const decltype(t1_ix_iter) iter1, const decltype(t2_ix_iter) iter2) {
      for (int i = 0; i != Ni; ++i)
        if (iter1[i] != iter2[i])
          return (iter1[i] < iter2[i]);
      return false;
    };
                      
  while (t1_ix_iter != t1_ix_iter_end &&
         t2_ix_iter != t2_ix_iter_end) {

    if (ixless(t1_ix_iter, t2_ix_iter)) {
      // multiplication by 0 -> no new coefficient, but advance index
      t1_ix_iter += Ni;
      ++t1_coef_iter;
    } else if (ixless(t2_ix_iter, t1_ix_iter)) {
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
  
  auto t1_ix_iter = ixs1.begin();
  auto t1_ix_iter_end = ixs1.end();
  auto t1_coef_iter = t1.coefs().begin();

  auto t2_ix_iter = ixs2.begin();
  auto t2_ix_iter_end = ixs2.end();
  auto t2_coef_iter = t2.coefs().begin();

  std::vector<T> target_coefs;
  std::vector<typename TensorComp<T>::Index> target_ixs;

  const auto ixless =
    [Ni](const decltype(t1_ix_iter) iter1, const decltype(t2_ix_iter) iter2) {
      for (int i = 0; i != Ni; ++i)
        if (iter1[i] != iter2[i])
          return (iter1[i] < iter2[i]);
      return false;
    };
                      
  while (t1_ix_iter != t1_ix_iter_end &&
         t2_ix_iter != t2_ix_iter_end) {

    if (ixless(t1_ix_iter, t2_ix_iter)) {
      // division by zero
      throw std::runtime_error("Division by zero in element-wise tensor division.");
      
    } else if (ixless(t2_ix_iter, t1_ix_iter)) {
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
