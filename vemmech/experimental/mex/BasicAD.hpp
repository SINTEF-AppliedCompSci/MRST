#ifndef __BASIC_AD_HPP
#define __BASIC_AD_HPP

#include <algorithm>
#include <vector>
#include <iostream>
// #include <string> // @@debug
#include <assert.h>

class BasicAD
{
public:
  BasicAD() : val_(0), num_derivs_(0) {}; 
  BasicAD(double val) :val_(val), num_derivs_(0) {}; 
  BasicAD(double val, int num_derivs) : val_(val), num_derivs_(num_derivs) {}

  BasicAD(const BasicAD& rhs) : val_(rhs.val()),
                                deriv_vals_(rhs.deriv_vals_),
                                deriv_ixs_(rhs.deriv_ixs_),
                                num_derivs_(rhs.num_derivs_) {}

  void swap(BasicAD& rhs) {std::swap(val_, rhs.val_);
                           std::swap(deriv_vals_, rhs.deriv_vals_);
                           std::swap(deriv_ixs_, rhs.deriv_ixs_);
                           std::swap(num_derivs_, rhs.num_derivs_);}
  
  void operator=(double val) {val_ = val;
                              deriv_vals_.clear();
                              deriv_ixs_.clear();
                              num_derivs_ = 0;}
  
  void operator=(const BasicAD& rhs) {BasicAD tmp(rhs); this->swap(tmp);}

  double val() const { return val_;}
  
  const std::vector<double> derivs() const {
    std::vector<double> derivs(num_derivs_, 0);
    assert(deriv_vals_.size() == deriv_ixs_.size());
    for (size_t i = 0; i != deriv_ixs_.size(); ++i) {
      assert(deriv_ixs_[i] < derivs.size());
      derivs[deriv_ixs_[i]] = deriv_vals_[i];
    }
    return derivs;
  }

  void setDeriv(size_t ix, double value) {
    if (ix >= num_derivs_)
      // increase number of derivatives
      num_derivs_ = ix+1;

    // later indices will overwrite earlier indices when calling 'derivs()' later
    deriv_ixs_.push_back(ix);
    deriv_vals_.push_back(value);
  }
      
  // addition/subtraction/multiplication/division with scalar
  void operator+=(double val) {val_ += val; }
  void operator-=(double val) {val_ -= val; }
  void operator*=(double val) {val_ *= val; mul_derivs_(val);}
  void operator/=(double val) {val_ /= val; mul_derivs_(1.0/val);}
  
  // addition/subtraction with BasicAD
  void operator+=(const BasicAD& rhs) { val_ += rhs.val(); add_derivs_(rhs); }
  void operator-=(const BasicAD& rhs) { val_ -= rhs.val(); add_derivs_(rhs, -1); }

  // multiplication/division
  void operator*=(const BasicAD& rhs) {
    mul_derivs_(rhs.val());
    add_derivs_(rhs, val_);
    val_ *= rhs.val();
    //dump_derivs("*=AD");
  }

  void operator/=(const BasicAD& rhs) {
    mul_derivs_(rhs.val());
    add_derivs_(rhs, -val_);
    mul_derivs_(1/(rhs.val() * rhs.val()));
    val_ /= rhs.val();
    //dump_derivs("/=AD");
  }

  void write(std::ostream& os) const {
    os << val_ << '\n';
    const auto ders = derivs();
    for (size_t i = 0; i != ders.size(); ++i)
      os << ders[i] << " ";
    os << '\n';
  }
  
private:
  double val_;
  std::vector<double> deriv_vals_;
  std::vector<size_t> deriv_ixs_;
  size_t num_derivs_;

  // inline void dump_derivs(std::string message) {
  //   std::cout << message << " - Value is: " << val_ << std::endl;
  //   for (int i = 0; i != derivs_.size(); ++i) 
  //     std::cout << derivs_[i] << " ";
  //   std::cout << std::endl;
  // }
  
  inline void mul_derivs_(double val) {
    std::for_each(deriv_vals_.begin(), deriv_vals_.end(), [val] (double& d) {d *= val;});
  }

  inline void add_derivs_(const BasicAD& rhs, const double factor=1) {

    std::vector<double> derivs = this->derivs();
    std::vector<double> other_derivs = rhs.derivs();
    if (derivs.size() < other_derivs.size())
      derivs.insert(derivs.end(), other_derivs.size() - derivs.size(), 0);
    assert(derivs.size() >= other_derivs.size());
    
    auto target = derivs.begin();
    for (auto it = other_derivs.begin(); it != other_derivs.end(); ++it, ++target)
      *target += *it * factor;

    // resetting internal data related to derivatives
    num_derivs_ = (num_derivs_ > rhs.num_derivs_) ? num_derivs_ : rhs.num_derivs_;
    assert(num_derivs_ == derivs.size());
    deriv_vals_.clear();
    deriv_ixs_.clear();
    for (size_t i = 0; i != derivs.size(); ++i) {
      if (derivs[i] != 0) {
        deriv_vals_.push_back(derivs[i]);
        deriv_ixs_.push_back(i);
      }
    }
  }
};

inline std::ostream& operator<<(std::ostream& os, const BasicAD& ad) {
  ad.write(os);
  return os;
}

inline BasicAD operator+(const BasicAD& a1, const BasicAD& a2) {
  BasicAD result(a1);
  result += a2;
  return result;
}

inline BasicAD operator-(const BasicAD& a1, const BasicAD& a2) {
  BasicAD result(a1);
  result -= a2;
  return result;
}

inline BasicAD operator*(const BasicAD& a1, const BasicAD& a2) {
  BasicAD result(a1);
  result *= a2;
  return result;
}

inline BasicAD operator/(const BasicAD& a1, const BasicAD& a2) {
  BasicAD result(a1);
  result /= a2;
  return result;
}

// template<typename T>
// inline BasicAD& operator+(const T& a1, const BasicAD& a2) {
//   BasicAD result(a1);
//   result += a2;
//   return result;
// }

// template<typename T>
// inline BasicAD& operator-(const T& a1, const BasicAD& a2) {
//   BasicAD result(a1);
//   result -= a2;
//   return result;
// }

// template<typename T>
// inline BasicAD& operator*(const T& a1, const BasicAD& a2) {
//   BasicAD result(a1);
//   result *= a2;
//   return result;
// }

// template<typename T>
// inline BasicAD& operator/(const T& a1, const BasicAD& a2) {
//   BasicAD result(a1);
//   result /= a2;
//   return result;
//}

#endif
