#pragma once
#include <functional>
#include <map>
#include <list>
#include <utility>
#include <numeric>

#include <cassert>
#include <iostream>

namespace kaspar {

class FixedHistogram {
public:
  FixedHistogram(int _B, std::vector<double> const& points, double eps=1e-5);
  
  double sum(double const& b); 

  void update(double const& p, int c=1); 

  void merge(FixedHistogram const& h);
  
  std::vector<double> uniform(int _B);
private:
  int B;
  double eps;
  std::map<double, int> data;
  double min_point, max_point;

  void _reduce(void); 

};

} // kaspar
