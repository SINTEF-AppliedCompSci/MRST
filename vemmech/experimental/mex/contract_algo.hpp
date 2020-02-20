#ifndef __CONTRACT_ALGO_HPP
#define __CONTRACT_ALGO_HPP

#include <set>
#include <numeric>
#include <algorithm>
#include <array>
#include <assert.h>

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
inline void advance_multiindex(std::vector<int>& mix, const std::vector<size_t>& bounds)
// ----------------------------------------------------------------------------
{
  mix.back() = mix.back() + 1;
  for (int i = mix.size() - 1; i > 0; --i) {
    if (mix[i] == bounds[i]) {
      mix[i] = 0;
      mix[i-1] = mix[i-1] + 1;
    }
  }
  std::cout << "new mix: ";
  for (int i = 0; i != mix.size(); ++i)
    std::cout << mix[i] << "  (" << bounds[i] << ") ";
  std::cout <<  std::endl;
}
  

// ----------------------------------------------------------------------------
template<typename T>
std::vector<std::vector<typename TensorComp<T>::Index>>
compute_multiindices(
     const std::vector<std::string>& cixs,
     const std::vector<int>& mix_comp,
     const std::vector<TensorComp<T>>& comps,
     std::vector<std::vector<std::array<typename TensorComp<T>::Index, 2>>>& ranges)
// ----------------------------------------------------------------------------
{
  typedef typename TensorComp<T>::Index Index;
  typedef std::array<Index, 2> RangeEntry;

  struct IxControl {
    int comp_ix;
    std::vector<Index> ix_entries;
    std::vector<std::vector<Index>> ix_values;
  };

  std::vector<IxControl> controls;
  
  for (int c_ix = 0; c_ix != comps.size(); ++c_ix) {

    // determine indices for which this component will control the multiindex
    std::vector<std::string> ixnames;
    std::vector<Index> ix_pos;
    for (int i = 0; i != mix_comp.size(); ++i) {
      if (mix_comp[i] == c_ix) {
        ixnames.push_back(cixs[i]);
        ix_pos.push_back(i);
      }
    }
    if (ixnames.empty())
      continue;
    
    std::vector<std::vector<Index>> ixvals(ixnames.size());
    for (int i = 0; i != ixnames.size(); ++i) {
      ixvals[i] = comps[c_ix].indexValuesFor(ixnames[i]);
      std::cout << "ixvals size: ";
      for (int j = 0; j != ixvals[i].size(); ++j)
        std::cout << ixvals[i][j] <<  " " ;
      std::cout << std::endl;
    }
    
    // identify where index changes occur
    std::vector<Index> changes(1, 0);
    for (int i = 0; i != ixvals[0].size(); ++i)
      for (int j = 0; j != ixvals.size() && changes.back() != i; ++j)
        if (ixvals[j][i-1] != ixvals[j][i]) 
          changes.push_back(i);

    std:: cout << "changes size: " << changes.size() << std::endl;
    
    std::vector<std::vector<Index>> unique_ixvals;
    for (int i = 0; i != ix_pos.size(); ++i) {
      std::vector<Index> unique_ival(changes.size());

      for (int j = 0; j != changes.size(); ++j)
        unique_ival[j] = ixvals[i][changes[j]];

      unique_ixvals.push_back(unique_ival);
    }

    controls.emplace_back(IxControl {c_ix, ix_pos, unique_ixvals} );
    std:: cout << "unique ixvals: " << unique_ixvals[0].size() << std::endl;
  }

  // sort controls so that highest stride lies in first index
  std::sort(controls.begin(), controls.end(),
            [](const IxControl& c1, const IxControl& c2) {
              return c1.ix_entries[0] < c2.ix_entries[0];});
  
  std::vector<int> indices(controls.size(), 0); 
  std::vector<size_t> indices_max;
  for (const auto& c : controls)
    indices_max.push_back(c.ix_values[0].size());
  assert(indices.size() > 0);
  
  std::vector<std::vector<Index>> result;
  std::vector<Index> new_multiindex(mix_comp.size());
  
  while (indices[0] != indices_max[0]) {

    for (int i = 0; i != controls.size(); ++i) {
      const auto& ctrl = controls[i];
      for (int j = 0; j != ctrl.ix_values.size(); ++j)
        new_multiindex[ctrl.ix_entries[j]] = ctrl.ix_values[j][indices[i]];
    }
    result.push_back(new_multiindex);
    
    // advance loop multiindex
    advance_multiindex(indices, indices_max);
  }

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
  std::vector<std::vector<std::array<Index, 2>>> comp_ranges;
  const std::vector<std::vector<Index>> mix =
    compute_multiindices(cixs, mix_comp, comps, comp_ranges);

  cout << "Start of multiindex list: " << endl;
  for (int i = 0; i != mix.size(); ++i) {
    const auto& entry = mix[i];
    for (int j = 0; j != mix[i].size(); ++j)
      cout << mix[i][j] << " ";
    cout << endl;
  }
  
  return comps; // @@ implement properly
}
  
#endif
