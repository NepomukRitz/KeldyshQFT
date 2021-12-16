# ifndef  DOMAIN_HPP
# define DOMAIN_HPP

#pragma once

#include <assert.h>

/*
template <typename T>
class Domain {
 private:
  virtual std::vector<Domain> split();
  virtual Base<T> size();
  virtual ~Domain();
};
*/

template <typename Q, typename Integrand>
class Domain1D {
 public:
  // double left_;
  // double right_;
  typedef Q value_type ;
  typedef Integrand f_type;
  //using value_type = T;
  //using f_type = F; //std::function<T(Base<T>)>;

  Domain1D(double left, double right) : left_(left), right_(right) {
    assert(left <= right);
  }

  std::vector<Domain1D> split() const {
    auto length = right_ - left_;
    auto middle = left_ + length / 2;

    return std::vector<Domain1D>({{left_, middle}, {middle, right_}});
  }

  template<typename TransformableIntegrand>
  auto transform(TransformableIntegrand &f) const {
    return  [&](double x) {
      double dmc = 0.5 * (right_ - left_);
      double dpc = 0.5 * (right_ + left_);

      Q val = f(x * dmc + dpc);
      //std::cout << "|" << x * dmc + dpc << "|" << val << "|\n";

      return val * dmc;
    };
  }

  double size() const { return right_ - left_; }

 private:
    double left_;
    double right_;

};

#endif // DOMAIN_HPP