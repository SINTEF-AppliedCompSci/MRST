#ifndef __BASIC_AD_HPP
#define __BASIC_AD_HPP

#include <algorithm>
#include <vector>
// #include <iostream> //@@ debug
// #include <string> // @@debug

class BasicAD
{
public:
  BasicAD() : val_(0) {};
  BasicAD(double val) :val_(val) {};
  BasicAD(double val, int num_derivs) : val_(val), derivs_(num_derivs, 0) {}

  BasicAD(const BasicAD& rhs) : val_(rhs.val()), derivs_(rhs.derivs()) {}

  void swap(BasicAD& rhs)            {std::swap(val_, rhs.val_); std::swap(derivs_, rhs.derivs_);}
  void operator=(double val)         {val_ = val; derivs_.clear();}
  void operator=(const BasicAD& rhs) {BasicAD tmp(rhs); this->swap(tmp);}

  double val() const { return val_;}
  const std::vector<double>& derivs() const { return derivs_;}

  void setDeriv(size_t ix, double value) {
    if (ix >= derivs_.size())
      derivs_.insert(derivs_.end(), ix-derivs_.size() + 1, 0);
    derivs_[ix] = value;
  }
      
  // addition/subtraction/multiplication/division with scalar
  void operator+=(double val) {val_ += val; }
  void operator-=(double val) {val_ -= val; }
  void operator*=(double val) {
    val_ *= val; std::for_each(derivs_.begin(), derivs_.end(), [val](double& d) {d *= val;}); }
  void operator/=(double val) {
    val_ /= val; std::for_each(derivs_.begin(), derivs_.end(), [val](double& d) {d /= val;}); }
  
  // addition/subtraction with BasicAD
  void operator+=(const BasicAD& rhs) { val_ += rhs.val(); add_derivs_(rhs.derivs()); }
  void operator-=(const BasicAD& rhs) { val_ -= rhs.val(); add_derivs_(rhs.derivs(), -1); }

  // multiplication/division
  void operator*=(const BasicAD& rhs) {
    mul_derivs_(rhs.val());
    add_derivs_(rhs.derivs(), val_);
    val_ *= rhs.val();
    //dump_derivs("*=AD");
  }

  void operator/=(const BasicAD& rhs) {
    mul_derivs_(rhs.val());
    add_derivs_(rhs.derivs_, -val_);
    mul_derivs_(1/(rhs.val() * rhs.val()));
    val_ /= rhs.val();
    //dump_derivs("/=AD");
  }
  
  
private:
  double val_;
  std::vector<double> derivs_;

  // inline void dump_derivs(std::string message) {
  //   std::cout << message << " - Value is: " << val_ << std::endl;
  //   for (int i = 0; i != derivs_.size(); ++i) 
  //     std::cout << derivs_[i] << " ";
  //   std::cout << std::endl;
  // }
  
  inline void mul_derivs_(double val) {
    std::for_each(derivs_.begin(), derivs_.end(), [val] (double& d) {d *= val;});
  }

  inline void add_derivs_(const std::vector<double>& others, const double factor=1) {
    if (derivs_.size() < others.size())
      derivs_.insert(derivs_.end(), others.size() - derivs_.size(), 0);

    auto target = derivs_.begin();
    for (auto it = others.begin(); it != others.end(); ++it, ++target)
      *target += *it * factor;
  }
};

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
