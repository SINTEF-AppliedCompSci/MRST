#ifndef __MULTIIX_TOOLS_HPP
#define __MULTIIX_TOOLS_HPP

#include <set>
#include <numeric>
#include <algorithm>
#include <array>
#include <assert.h>

// ----------------------------------------------------------------------------
template<typename Index>
std::vector<size_t> sort_multiindex(std::vector<Index>& mix, int num_elem)
// ----------------------------------------------------------------------------
{
  struct SortElem {
    const Index* ptr;
    size_t init_pos;
  };
  size_t N = mix.size() / num_elem;
  assert(num_elem * N == mix.size());
  
  std::vector<SortElem> mixptrs(N);
  for (size_t i = 0; i != N; ++i)
    mixptrs[i] = { &mix[i * num_elem], i};

  auto comp = [num_elem](const SortElem& i1, const SortElem& i2) {
                for (int i = 0; i != num_elem; ++i)
                  if (i1.ptr[i] != i2.ptr[i])
                    return i1.ptr[i] < i2.ptr[i];
                return false;
              };
  std::sort(mixptrs.begin(), mixptrs.end(), comp);

  // copying sorted elements back
  std::vector<size_t> reindexes(N);
  std::vector<Index> result(mix.size());
  size_t count = 0;
  size_t reindex_count = 0;
  for (auto it = mixptrs.begin(); it != mixptrs.end(); ++it) {
    for (int i = 0; i != num_elem; ++i) 
      result[count++] = it->ptr[i];
    
    reindexes[reindex_count++] = it->init_pos;
  }
  assert(count == result.size());
  assert(reindex_count == reindexes.size());
  
  std::swap(mix, result);
  return reindexes;
}

// ----------------------------------------------------------------------------
inline void advance_multiindex(std::vector<size_t>& mix, const std::vector<size_t>& bounds)
// ----------------------------------------------------------------------------
{
  assert(mix.size() > 0);
  assert(bounds.size() == mix.size());
  
  mix.back() = mix.back() + 1;
  for (int i = mix.size() - 1; i > 0; --i) {
    if (mix[i] == bounds[i]) {
      mix[i] = 0;
      mix[i-1] = mix[i-1] + 1;
    }
  }
}


// ----------------------------------------------------------------------------
template<typename Index>
void transpose(std::vector<Index>& vals, size_t num_cols)
// ----------------------------------------------------------------------------
{
  if (num_cols == 0)
    // nothing to do
    return;
  
  const size_t num_rows = vals.size() / num_cols;
  assert(num_rows * num_cols == vals.size());
  std::vector<Index> result(vals.size());
  size_t ix = 0;
  for(size_t i = 0; i != num_rows; ++i)
    for(size_t j = 0; j != num_cols; ++j)
      result[ix++] = vals[j * num_rows + i];

  assert(ix == result.size());
  std::swap(vals, result);
}

#endif
