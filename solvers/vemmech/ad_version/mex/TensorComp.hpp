#ifndef __TENSOR_COMP_HPP
#define __TENSOR_COMP_HPP

#include <iostream>
#include <vector>
#include <string>
#include <tuple>
#include <algorithm>
//#include <parallel/algorithm>  //@@ gcc_specific
#include <array>
#include <set>
#include <cmath>
#include <assert.h>

// ============================================================================
template<typename T> class TensorComp
// ============================================================================
{
public:

  typedef size_t Index;

  // empty tensor
  TensorComp() {}
  // intrinsic scalar
  TensorComp(const T& scalar) : indexnames_(std::vector<std::string>()),
                                coefs_(std::vector<T>(1, scalar)),
                                ixs_(std::vector<Index>()) {}
  // indexed tensor
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
  //TensorComp<T>& sortElementsByIndex_old(bool descending=false);
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
  std::string tmp;
  for (int i = 0; i != num; ++i) {
    is >> tmp;
    if (tmp == "-nan")
      coefs[i] = std::nan("");
    else
      coefs[i] = std::atof(tmp.c_str());
  }

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
    if (ixname == indexNames()[i]) {
      assert( (i+1) * numCoefs() <= ixs_.size());
      const auto startpos = ixs_.begin() + i * numCoefs();
      const auto endpos = ixs_.begin() + (i+1) * numCoefs();
      return std::vector<Index>(startpos, endpos);
    }      
  return std::vector<Index>();
}

// ----------------------------------------------------------------------------
template<typename T> inline size_t
TensorComp<T>::numUniqueValuesFor(const std::string& ixname) const
// ----------------------------------------------------------------------------
{
  const auto it = std::find(indexNames().begin(), indexNames().end(), ixname);

  return(it == indexNames().end()) ?
    0 :
    numUniqueValuesFor(it - indexNames().begin());
}

// ----------------------------------------------------------------------------
template<typename T> inline size_t
TensorComp<T>::numUniqueValuesFor(int ix) const
// ----------------------------------------------------------------------------
{
  const size_t N = numCoefs();

  assert(ix < numIndices());
  const auto startpos = ixs_.begin() + (ix*N);
  const auto endpos = ixs_.begin() + ((ix+1)*N);

  return std::set<Index>(startpos, endpos).size();
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
  //std::cout << "Entering sumEqualIndices. " << std::endl;
  sortElementsByIndex(); // ensure equal indices follow each other

  //std::cout << "Still in sumEqualIndices - just finished sorting elements." << std::endl;
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

  // @@ THIS LOOP SHOULD BE PARALLELIZED
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

  //std::cout << "Finishing sumEqualIndices." << std::endl;
  return *this;
}


// ----------------------------------------------------------------------------
template<typename T> inline TensorComp<T>&
TensorComp<T>::sortIndicesByNumber(bool descending)
// ----------------------------------------------------------------------------
{
  // index name included to ensure unique sort
  typedef std::tuple<size_t, std::string, size_t> ICount; 

  // determine unique entries in for each index
  std::vector<ICount> num_uniques(numIndices());

  const size_t N = numCoefs();
  for (int i = 0; i != numIndices(); ++i)
    num_uniques[i] = ICount {numUniqueValuesFor(i), indexnames_[i], i};

  // permute indices accordingly
  auto comp = descending ? [](const ICount& a, const ICount& b) { return a > b;} :
                           [](const ICount& a, const ICount& b) { return a < b;};

    // [](const ICount& a, const ICount& b) { return std::get<0>(a) > std::get<0>(b);} :
    // [](const ICount& a, const ICount& b) { return std::get<0>(a) < std::get<0>(b);};
  
  std::stable_sort(num_uniques.begin(), num_uniques.end(), comp);
  // if (descending)
  //   std::reverse(num_uniques.begin(), num_uniques.end());

  std::vector<int> perm(numIndices());
  for (int i = 0; i != numIndices(); ++i)
    perm[i] = std::get<2>(num_uniques[i]);

  return permuteIndices(perm);
}

// ----------------------------------------------------------------------------
template<typename T> inline TensorComp<T>&
TensorComp<T>::sortElementsByIndex(bool descending)
// ----------------------------------------------------------------------------
{
  const size_t Nc = numCoefs();
  const size_t Ni = numIndices();

  //std::cout << "Sorting " << Nc << " indices." << std::endl;
  std::vector<const Index*> ptrs(Nc); // first sort pointers, then shuffle data

  // make pointers point to first index element of each multiindex
  for (size_t i = 0; i != Nc; ++i)
    ptrs[i] = &ixs_[i];

  auto comp = [Ni, Nc] (const Index* i1, const Index* i2) {
                for (int i = 0; i != Ni; ++i, i1 += Nc, i2 += Nc)
                  if (*i1 != *i2)
                    return *i1 < *i2;
                return false;
              };

  if (std::is_sorted(ptrs.begin(), ptrs.end(), comp)) {
    //std::cout << "Already sorted. Returning." << std::endl;
    return *this;
  }
  
  std::sort(ptrs.begin(), ptrs.end(), comp);
  //__gnu_parallel::sort(ptrs.begin(), ptrs.end(), comp);

  if (descending)
    std::reverse(ptrs.begin(), ptrs.end());
  
  // reshuffle data (indices and coefficients)
  std::vector<T> coefs_sorted(coefs_.size());
  for (size_t ic = 0; ic != Nc; ++ic) {
    coefs_sorted[ic] = coefs_[size_t(ptrs[ic] - &ixs_[0])];
  }
  
  std::vector<Index> ixs_sorted(ixs_.size());
  auto iter = ixs_sorted.begin();
  for (size_t ii = 0; ii != Ni; ++ii) {
    for (size_t ic = 0; ic != Nc; ++ic) {
      *iter++ = *ptrs[ic];
      ptrs[ic] += Nc; // progress pointer to next column of indices
    }
  }
  //std::cout << "Finished!" << std::endl;
  std::swap(coefs_, coefs_sorted);
  std::swap(ixs_, ixs_sorted);
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
  std::vector<Index> result(ixs_.size(), Index(0));
  std::vector<std::string> ixnames_new(numIndices());
  
  const size_t N = numCoefs();
  for (int i = 0; i != numIndices(); ++i) {
    const auto startpos = ixs_.begin() + perm[i] * N;
    const auto endpos = ixs_.begin() + (perm[i] + 1) * N;
    std::copy(startpos, endpos, &result[i*N]);
    ixnames_new[i] = indexNames()[perm[i]];
  }
  std::swap(ixs_, result);
  std::swap(indexnames_, ixnames_new);

  return *this;
}
  
#endif
