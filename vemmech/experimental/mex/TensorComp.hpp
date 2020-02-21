#ifndef __TENSOR_COMP_HPP
#define __TENSOR_COMP_HPP

#include <iostream>
#include <vector>
#include <string>
#include <tuple>
#include <algorithm>
#include <array>
#include <set>
#include <assert.h>

// ============================================================================
template<typename T> class TensorComp
// ============================================================================
{
public:
  //using Index = int;
  using ICount = std::tuple<int, size_t>;

  typedef int Index;
  TensorComp() {}
  TensorComp(const std::vector<std::string>& indexnames,
             const std::vector<T>& coefs,
             const std::vector<Index>& ixs) :
    indexnames_(indexnames), coefs_(coefs), ixs_(ixs) {}

  const std::vector<std::string>& indexNames() const {return indexnames_;}
  const size_t numIndices() const {return indexnames_.size();}
  const size_t numCoefs() const {return coefs_.size();}
  const std::vector<T>& coefs() const {return coefs_;}
  const std::vector<Index>& ixs() const { return ixs_;}

  // return position of index with given name, or -1 if it doesn't exist
  int indexPos(const std::string& name) const;

  std::vector<Index> indexValuesFor(const std::string& ixname) const;
  
  template<typename Indexable>
  TensorComp<T>& permuteIndices(const Indexable perm);
  TensorComp<T>& sortElementsByIndex(bool descending=false);
  TensorComp<T>& sortIndicesByNumber(bool descending=false);
  TensorComp<T>& sumEqualIndices();
  
  
  size_t numUniqueValuesFor(int ix) const ;
  size_t numUniqueValuesFor(const std::string& ixname) const;
  
  // any index name found in 'ixnames' should be moved to the front of the
  // multiindex
  TensorComp<T>& moveIndicesFirst(const std::vector<std::string>& ixnames);

  void write(std::ostream& os) const ;
  void read(std::istream& is);
  
private:
  std::vector<std::string> indexnames_;
  std::vector<T> coefs_;
  std::vector<Index> ixs_;

  template<typename Indexable>
  static bool is_permutation_(const Indexable perm, const Index num);
  
}; // end class TensorComp

// ----------------------------------------------------------------------------
template<typename T>
void TensorComp<T>::write(std::ostream& os) const
// ----------------------------------------------------------------------------
{
  // writing index names
  os << indexnames_.size() << '\n';
  for (const auto& n : indexnames_)
    os << n << '\n';

  // writing coefficients
  os << coefs_.size() << '\n';
  for (const auto& c : coefs_)
    os << c << ' ';
  os << '\n';

  // writing indices
  os << ixs_.size() << '\n';
  for (const auto& i : ixs_)
    os << i << ' ';
  os << '\n';
}

// ----------------------------------------------------------------------------
template<typename T>
void TensorComp<T>::read(std::istream& is)
// ----------------------------------------------------------------------------
{
  // read index names
  int num; is >> num;
  std::vector<std::string> ixnames(num);
  
  for (int i = 0; i != num; ++i)
    is >> ixnames[i];

  // read coefficients
  is >> num;
  std::vector<T> coefs(num);
  for (int i = 0; i != num; ++i)
    is >> coefs[i];

  // read indices
  is >> num;
  std::vector<Index> ixs(num);
  for (int i = 0; i != num; ++i)
    is >> ixs[i];

  indexnames_.swap(ixnames);
  coefs_.swap(coefs);
  ixs_.swap(ixs);
}


// ----------------------------------------------------------------------------
template<typename T> inline 
int TensorComp<T>::indexPos(const std::string& name) const
// ----------------------------------------------------------------------------
{
  const auto pos = std::find(indexnames_.begin(), indexnames_.end(), name);
  return (pos == indexnames_.end()) ? -1 : pos - indexnames_.begin();
}
  
// ----------------------------------------------------------------------------
template<typename T> inline std::vector<typename TensorComp<T>::Index>
TensorComp<T>::indexValuesFor(const std::string& ixname) const
// ----------------------------------------------------------------------------
{
  for (int i = 0; i != indexNames().size(); ++i)
    if (ixname == indexNames()[i])
      return std::vector<Index>(&ixs_[i * numCoefs()], &ixs_[(i+1) * numCoefs()]);
  return std::vector<Index>();
}

// ----------------------------------------------------------------------------
template<typename T> inline size_t
TensorComp<T>::numUniqueValuesFor(const std::string& ixname) const
// ----------------------------------------------------------------------------
{
  const auto it = std::find(indexNames().begin(), indexNames().end(), ixname);

  return (it == indexNames().end()) ?
    0 :
    numUniqueValuesFor(it - indexNames().begin());
}

// ----------------------------------------------------------------------------
template<typename T> inline size_t
TensorComp<T>::numUniqueValuesFor(int ix) const
// ----------------------------------------------------------------------------
{
  const size_t N = numCoefs();
  return std::set<Index>(&ixs_[ix*N], &ixs_[(ix+1)*N]).size();
}

// ----------------------------------------------------------------------------
template<typename T> inline TensorComp<T>&
TensorComp<T>::moveIndicesFirst(const std::vector<std::string>& ixnames)
// ----------------------------------------------------------------------------
{
  std::vector<int> perm, keep_ix;
  for (int i = 0; i != numIndices(); ++i)
    if (std::find(ixnames.begin(), ixnames.end(), indexNames()[i]) != ixnames.end())
      perm.push_back(i);
    else
      keep_ix.push_back(i);

  perm.insert(perm.end(), keep_ix.begin(), keep_ix.end());

  return permuteIndices(perm);
  
}

// ----------------------------------------------------------------------------
template<typename T> inline TensorComp<T>& TensorComp<T>::sumEqualIndices()
// ----------------------------------------------------------------------------
{
  sortElementsByIndex(); // ensure equal indices follow each other

  const int N = numIndices();
  std::vector<T> new_coefs;
  std::vector<std::vector<Index>> new_indices(N);

  auto coef_ptr = coefs().begin();
  std::vector<const Index*> old_ix_ptrs(N);

  new_coefs.push_back(*coef_ptr++);
  for (int i = 0; i != N; ++i) {
    old_ix_ptrs[i] = &(ixs_[i * coefs().size()]);
    new_indices[i].push_back(*old_ix_ptrs[i]++);
  }
  
  while(coef_ptr != coefs().end()) {
    // check if current index set equals 
    bool same = true;
    for (int i = 0; i != N && same; ++i)
      same = *(old_ix_ptrs[i]) == new_indices[i].back();

    if (same) {
      new_coefs.back() += *coef_ptr; // add to existing element
    } else {
      new_coefs.push_back(*coef_ptr);
      for (int i = 0; i != N; ++i)
        new_indices[i].push_back(*old_ix_ptrs[i]);
    }

    // advance iterators
    ++coef_ptr;
    for (int i = 0; i != N; ++i)
      ++old_ix_ptrs[i];
  }

  std::vector<Index> all_ixs;
  for (int i = 0; i != N; ++i)
    all_ixs.insert(all_ixs.end(), new_indices[i].begin(), new_indices[i].end());

  // replace old member data with new
  coefs_.swap(new_coefs);
  ixs_.swap(all_ixs);

  return *this;
}


// ----------------------------------------------------------------------------
template<typename T> inline TensorComp<T>&
TensorComp<T>::sortIndicesByNumber(bool descending)
// ----------------------------------------------------------------------------
{
  // determine unique entries in for each index
  std::vector<ICount> num_uniques(numIndices());

  const size_t N = numCoefs();
  for (int i = 0; i != numIndices(); ++i)
    num_uniques[i] = ICount {numUniqueValuesFor(i), i};

  // permute indices accordingly
  std::sort(num_uniques.begin(), num_uniques.end());
  if (descending)
    std::reverse(num_uniques.begin(), num_uniques.end());

  std::vector<int> perm(numIndices());
  for (int i = 0; i != numIndices(); ++i)
    perm[i] = std::get<1>(num_uniques[i]);

  return permuteIndices(perm);
  
}


// ----------------------------------------------------------------------------
template<typename T> inline TensorComp<T>&
TensorComp<T>::sortElementsByIndex(bool descending)
// ----------------------------------------------------------------------------
{
  using Ivec = std::vector<Index>;
  Ivec elem(numIndices()+1, 0);
  std::vector<Ivec> tmp_sort(numCoefs(), elem);
  const size_t N = numCoefs();
  for (int cix = 0; cix != numCoefs(); ++cix) {
    for (int i = 0; i != numIndices(); ++i)
      tmp_sort[cix][i] = ixs_[i * N + cix];
    tmp_sort[cix][numIndices()] = cix;
  }

  std::sort(tmp_sort.begin(), tmp_sort.end(),
            [](const Ivec& v1, const Ivec& v2) {
              for (int i = 0; i != v1.size(); ++i)
                if (v1[i] != v2[i])
                  return (v1[i] < v2[i]);
              return true;
            });

  if (descending)
    std::reverse(tmp_sort.begin(), tmp_sort.end());
  
  // copy sorted indices back
  for (int cix = 0; cix != numCoefs(); ++cix)
    for (int i = 0; i != numIndices(); ++i)
      ixs_[i * N + cix] = tmp_sort[cix][i];
  
  // sorting coefs
  std::vector<T> coefs_sorted(coefs_.size());
  for (int cix = 0; cix != numCoefs(); ++cix)
    coefs_sorted[cix] = coefs_[tmp_sort[cix][numIndices()]];

  std::swap(coefs_, coefs_sorted);

  return *this;
  
}
  
// ----------------------------------------------------------------------------
template<typename T> template<typename Indexable> inline bool
TensorComp<T>::is_permutation_(const Indexable perm, const Index num)
// ----------------------------------------------------------------------------
{
  std::vector<Index> input(&perm[0], &perm[0] + num);
  std::sort(input.begin(), input.end());
  for (int i = 0; i != input.size(); ++i)
    if (input[i] != i)
      return false;
  return true;
}

// ----------------------------------------------------------------------------
template<typename T> template<typename Indexable> inline TensorComp<T>&
TensorComp<T>::permuteIndices(const Indexable perm)
// ----------------------------------------------------------------------------
{
  assert(is_permutation_(perm, numIndices()));
  std::vector<Index> result(ixs_.size(), T(0));
  std::vector<std::string> ixnames_new(numIndices());
  
  const size_t N = numCoefs();
  for (int i = 0; i != numIndices(); ++i) {
    std::copy(&ixs_[perm[i] * N], &ixs_[(perm[i] + 1) * N], &result[i*N]);
    ixnames_new[i] = indexNames()[perm[i]];
  }
  std::swap(ixs_, result);
  std::swap(indexnames_, ixnames_new);

  return *this;
}
  
#endif
