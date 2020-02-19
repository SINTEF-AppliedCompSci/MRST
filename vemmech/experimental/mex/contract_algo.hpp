#ifndef __CONTRACT_ALGO_HPP
#define __CONTRACT_ALGO_HPP

#include <set>
#include <numeric>
#include <algorithm>

#include "TensorComp.hpp"

// ----------------------------------------------------------------------------
template<typename T> inline std::vector<std::string>
identify_contracting_multiindex(const std::vector<TensorComp<T>>& comps)
// ----------------------------------------------------------------------------
{
  // extract all indices
  std::vector<std::string> all_ixs;
  for (const auto& c : comps)
    std::copy(c.indexNames().begin(),
              c.indexNames().end(),
              std::back_inserter(all_ixs));

  // identify repeated indices
  std::set<std::string> contr_ixs;
  std::sort(all_ixs.begin(), all_ixs.end());
  for (int i = 0; i != all_ixs.size()-1; ++i)
    if (all_ixs[i] == all_ixs[i+1])
      contr_ixs.insert(all_ixs[i]);

  // converting into a vector
  const std::vector<std::string> cixs(contr_ixs.begin(), contr_ixs.end());

  
  // determine order of indices (start index should have highest number of
  // unique values)
  std::vector<size_t> ix_numvals(cixs.size(), 0);
  for (int i = 0; i != cixs.size(); ++i)
    for (const auto& c : comps)
      ix_numvals[i] = std::max(ix_numvals[i], c.numUniqueValuesFor(cixs[i]));

  std::vector<std::array<int, 2>> sort_ix;
  for (int i = 0; i != ix_numvals.size(); ++i)
    sort_ix.emplace_back( std::array<int, 2>{(int)ix_numvals[i], i} );

  std::sort(sort_ix.begin(), sort_ix.end());
  std::reverse(sort_ix.begin(), sort_ix.end());

  std::vector<int> perm(sort_ix.size());
  std::transform(sort_ix.begin(), sort_ix.end(), perm.begin(),
                 [](const std::array<int, 2>& el) {return el[1];});
    
  std::vector<std::string> sorted_ixs(cixs.size());
  for (int i = 0; i != cixs.size(); ++i)
    sorted_ixs[i] = cixs[perm[i]];
  //return cixs;
  return sorted_ixs;
  
}

// ----------------------------------------------------------------------------
template<typename T> std::vector<int>
identify_comps_for_multiix_iteration(std::vector<TensorComp<T>>& comps,
                                    const std::vector<std::string>& cixs)
// ----------------------------------------------------------------------------
{
  std::vector<int> result(cixs.size(), -1); // -1 flags that index not yet set
  std::vector<int> count(cixs.size(), -1);

  for (int c_ix = 0; c_ix != comps.size(); ++c_ix) {
    for (int i = 0; i != cixs.size(); ++i) {
      const size_t numvals = comps[c_ix].numUniqueValuesFor(cixs[i]);
      if (numvals > 0 && (count[i] < 0 || count[i] > numvals)) {
        result[i] = c_ix;
        count[i] = numvals;
      }
    }
  }
  return result;
}

// ----------------------------------------------------------------------------
template<typename T>
std::vector<typename TensorComp<T>::Index>
compute_multiindices(
     const std::vector<std::string>& cixs,
     const std::vector<int>& mix_comp,
     const std::vector<TensorComp<T>>& comps,
     std::vector<std::vector<typename TensorComp<T>::Index>>& ranges)
// ----------------------------------------------------------------------------
{
  for (int c_ix = 0; c_ix != comps.size(); ++c_ix) {

    // determine indices for which this component will control the multiindex
    std::vector<std::string> ixs;
    for (int i = 0; i != mix_comp.size(); ++i)
      if (mix_comp[i] == c_ix)
        ixs.push_back(cixs[i]);

    std::vector<vector<Index>> ixvals(ixs.size());
    for (int i = 0; i != ixs.size(); ++i)
      ixvals[i] = comps[c_ix].indexValuesFor(ixs[i]);

    // identify where index changes occur
    vector<TensorComp<T>::Index> changes(1, 0);
    for (int i = 1; i != comps[c_ix].numCoefs(); ++i) 
      for (int j = 1; j != ixvals.size() && changes.back() != i; ++j) 
        if (ixvals[j][i-1] != ixvals[j][i]) 
          changes.push_back(i);

    ranges.push_back(changes);

    // @@ add multiindex stuff here
  }

  
  std::vector<typename TensorComp<T>::Index> result(1000);
  std::iota(result.begin(), result.end(), 1);
  return result;
}

// ============================================================================
template<typename T>
std::vector<TensorComp<T>> contract_components(std::vector<TensorComp<T>> comps)
// ============================================================================
{
  using std::cout;
  using std::endl;
  typedef typename TensorComp<T>::Index Index;
  // identify contracting indices; sort them according to number of unique values
    const std::vector<std::string> cixs(identify_contracting_multiindex(comps));

  cout << "Contracting indices: " << endl;
  for (auto l : cixs)
    cout << l << endl;

  // Determine which components will control the iteration of each index in the
  // multiindex.  
  const auto mix_comp = identify_comps_for_multiix_iteration(comps, cixs);

  cout << "Mix comp: " << endl;
  for (int i = 0; i != mix_comp.size(); ++i)
    cout << mix_comp[i] << endl;
  
  // prepare components for iteration  by sorting index order and coefficient properly
  for (int i = 0; i != comps.size(); ++i)  
    comps[i].sortIndicesByNumber(true).moveIndicesFirst(cixs).sortElementsByIndex();
    
  // identify all possible multiindex values and corresponding entry ranges in
  // components
  std::vector<std::vector<Index>> comp_ranges;
  const std::vector<Index> mix =
    compute_multiindices(cixs, mix_comp, comps, comp_ranges);

  cout << "Start of multiindex list: " << endl;
  for (int i = 0; i != 20; ++i) {
    for (int j = 0; j != cixs.size(); ++j) {
      cout << mix[i * cixs.size() + j] << " ";
    }
    cout << endl;
  }
  
  //auto result = comps[0];
  //result.permuteIndices(std::array<int, 4>{0,2,1});
  //result.sortElementsByIndex(true);
  // result.sortIndicesByNumber(true);

  // std::vector<std::string> tull {"dim"};

  // result.moveIndicesFirst(tull);
  
  return comps; // @@ implement properly
}
  
#endif
