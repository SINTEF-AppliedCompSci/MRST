#ifndef __TENSOR_COMP_HPP
#define __TENSOR_COMP_HPP

#include <vector>
#include <string>

// ============================================================================
template<typename T> class TensorComp
// ============================================================================
{
public:
  TensorComp(const std::vector<std::string>& indexnames,
             const std::vector<T>& coefs,
             const std::vector<int>& ixs) :
    indexnames_(indexnames), coefs_(coefs), ixs_(ixs) {}

  const std::vector<std::string>& indexNames() const {return indexnames_;}
  const size_t numIndices() const {return indexnames_.size();}
  const std::vector<T>& coefs() const {return coefs_;}
  const std::vector<int>& ixs() const { return ixs_;}
  
  
private:
  std::vector<std::string> indexnames_;
  std::vector<T> coefs_;
  std::vector<int> ixs_;
  
}; // end class TensorComp

// ----------------------------------------------------------------------------
template<typename T>
TensorComp<T> contract_components(const std::vector<TensorComp<T>>& comps)
// ----------------------------------------------------------------------------
{
  return comps[0]; // @@ implement properly
}
  




#endif
