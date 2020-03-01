#ifndef __BASIC_AD_HPP
#define __BASIC_AD_HPP

#include <algorithm>
#include <vector>

class BasicAD
{
public:
  BasicAD(double val) val_(val) {};
  BasicAD(double val, int num_derivs) : val_(val), derivs_(num_derivs, 0) {}

  BasicAD(const BasicAD& rhs) val_(rhs.val()), derivs_(rhs.derivs()) {}

  void swap(BasicAD& rhs)            {std::swap(val_, rhs.val()); std::swap(derivs_, rhs.derivs());}
  void operator=(double val)         {val_ = val; derivs_.clear();}
  void operator=(const BasicAD& rhs) {BasicAD tmp(rhs); this->swap(rhs);}

  double val() const { return val;}
  const std::ector<double>& derivs() const { return derivs_;}
  
  // addition/subtraction/multiplication/division with scalar
  void operator+=(val) {val_ += val;}
  void operator-=(val) {val_ -= val;}
  void operator*=(val) {val_ *= val; std::for_each(derivs_.begin(), derivs_.end(), [val](double& d) {d *= val;});}
  void operator/=(val) {val_ /= val; std::for_each(derivs_.begin(), derivs_.end(), [val](double& d) {d /= val;});}
  
  // addition/subtraction with BasicAD
  void operator+=(const BasicAD& rhs) { val_ += rhs.val(); add_derivs_(rhs.derivs());  }
  void operator-=(const BasicAD& rhs) { val_ -= rhs.val(); add_derivs_(rhs.derivs(), -1);}

  // multiplication/division
  void operator*=(const BasicAD& rhs) {
    val_ *= rhs.val();
    mul_derivs_(rhs.val());
    add_derivs_(rhs.derivs(), val_);
  }

  void operator/=(const BasicAD& rhs) {
    val_ /= rhs.val();
    mul_derivs_(rhs.val());
    add_derivs_(rhs.derivs, -val_);
    mul_derivs_(1/(rhs.val() * rhs.val()));
  }
  
  
private:
  double val_;
  std::vector<double> derivs_;

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

template<typename T>
inline BasicAD& operator+(const BasicAD& a1, const T& a2) {
  BasicAD result(a1);
  result += a2;
  return result;
}

template<typename T>
inline BasicAD& operator-(const BasicAD& a1, const T& a2) {
  BasicAD result(a1);
  result -= a2;
  return result;
}

template<typename T>
inline BasicAD& operator*(const BasicAD& a1, const T& a2) {
  BasicAD result(a1);
  result *= a2;
  return result;
}

template<typename T>
inline BasicAD& operator/(const BasicAD& a1, const T& a2) {
  BasicAD result(a1);
  result /= a2;
  return result;
}

template<typename T>
inline BasicAD& operator+(const T& a1, const BasicAD& a2) {
  BasicAD result(a1);
  result += a2;
  return result;
}

template<typename T>
inline BasicAD& operator-(const T& a1, const BasicAD& a2) {
  BasicAD result(a1);
  result -= a2;
  return result;
}

template<typename T>
inline BasicAD& operator*(const T& a1, const BasicAD& a2) {
  BasicAD result(a1)
  result *= a2;
  return result;
}

template<typename T>
inline BasicAD& operator/(const T& a1, const BasicAD& a2) {
  BasicAD result(a1);
  result /= a2;
  return result;
}

#endif
