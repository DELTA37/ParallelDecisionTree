#include <ParallelDecisionTree/FixedHistogram/fixed_histogram.h>
#include <cmath>

namespace kaspar {

void FixedHistogram::_reduce(void) {
  size_t size = data.size();

  std::list<double> diff;
  std::transform(data.begin(), std::prev(data.end()), std::next(data.begin()), std::back_inserter(diff), 
    [](std::pair<double, int> const& x, std::pair<double, int> const& y) -> double {
      return y.first - x.first;
    }
  );

  for (size_t i = 0; i < size - B; ++i) {
    auto list_it = std::min_element(diff.begin(), diff.end());
    auto list_nit = std::next(list_it);
    size_t ind = std::distance(diff.begin(), list_it);

    auto min_it = std::next(data.begin(), ind);
    auto min_nit = std::next(min_it);
    auto min_nnit = std::next(min_nit);
    auto min_pit = std::prev(min_it);

    double point = (min_it->first * min_it->second + min_nit->first * min_nit->second) / (min_it->second + min_nit->second);
    int value = min_it->second + min_nit->second;

    if (min_pit != data.end()) {
      diff.insert(list_it, point - min_pit->first);
    }
    if (min_nnit != data.end()) {
      diff.insert(list_it, min_nnit->first - point);
    }

    diff.erase(list_it);
    if (list_nit != diff.end()) {
      diff.erase(list_nit);
    }

    data[point] = value;
    data.erase(min_it);
    data.erase(min_nit);
  }
}

FixedHistogram::FixedHistogram(int _B, std::vector<double> const& points, double _eps) : B(_B), eps(_eps) {
  min_point = std::numeric_limits<double>::max();
  max_point = std::numeric_limits<double>::min();
  assert(_B == points.size());
  for (auto const& p : points) {
    data[p] += 1;
    min_point = std::min(p, min_point);
    max_point = std::max(p, max_point);
  }
}

double FixedHistogram::sum(double const& b) {
  auto nit = data.upper_bound(b);
  auto it = (nit != data.begin()) ? std::prev(nit) : data.end();
  double s = 0;
  double p_i, p_ni;
  int m_i, m_ni;
  std::tie(p_i, m_i) = (it != data.end()) ? *it : std::make_tuple(std::min(b, min_point - eps), 0);
  std::tie(p_ni, m_ni) = (nit != data.end()) ? *nit : std::make_tuple(std::max(max_point + eps, b), 0);

  double m_b = m_i + (m_ni - m_i) * (b - p_i) / (p_ni - p_i);
  s += (m_i + m_b) * (b - p_i) / (p_ni - p_i) / 2.0;

  for (auto _it = data.begin(); _it != it; ++_it) {
    s += _it->second;
  }
  s += m_i / 2.0;
  return s;
}

std::vector<double> FixedHistogram::uniform(int _B) {
  std::vector<double> s;
  std::vector<double> u;
  u.reserve(_B);
  s.reserve(B);
  // m is actually a FixedHistogram.sum([-inf, p_i])
  double sumed = 0;
  for (auto it = data.begin(); it != data.end(); ++it) {
    s.emplace_back(sumed + (it->first / 2.0));
    sumed += it->first;
  }
  u.emplace_back(min_point); 
  for (size_t j = 1; j < _B - 1; ++j) {
    double  b        = (sumed * j) / _B;
    auto    it       = std::prev(std::upper_bound(s.begin(), s.end(), b));
    auto    data_it  = std::next(data.begin(), std::distance(s.begin(), it));
    auto    data_nit = std::next(data_it);

    auto [p_i, m_i]   = *data_it;
    auto [p_ni, m_ni] = (data_nit != data.end()) ? *data_nit : std::tuple{max_point, 0};

    double s_i = *it;
    double d   = b - s_i;
    u.emplace_back(
      (std::abs(m_ni - m_i) < eps) ? p_i : (p_i + (p_ni - p_i) * (-m_i + std::sqrt(m_i * m_i + 2.0 * (m_ni - m_i) * d)) / (m_ni - m_i))
    );
  }
  u.emplace_back(max_point); 
  return u;
}

void FixedHistogram::update(double const& p, int c) {
  data[p] += c;
  min_point = std::min(p, min_point);
  max_point = std::max(p, max_point);
  _reduce();
}

void FixedHistogram::merge(FixedHistogram const& h) {
  for (auto it = h.data.begin(); it != h.data.end(); ++it) {
    data[it->first] += it->second;
    min_point = std::min(it->first, min_point);
    max_point = std::max(it->first, max_point);
  }
  _reduce();
}


} // kaspar

