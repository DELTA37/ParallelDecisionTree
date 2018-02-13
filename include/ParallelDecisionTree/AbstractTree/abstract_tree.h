#pragma once
#include <type_traits>
#include <vector>
#include <algorithm>
#include <array>
#include <list>

namespace kaspar {

template <typename XType, typename YType, typename RuleType, int N>
class AbstractTree {
public:
  typedef std::pair<XType, YType> ExampleType;

protected:
  struct Node {
    virtual ~Node(void) {}
  };
  struct RuleNode : public Node {
    RuleType rule;
    std::array<typename std::list<Node*>::iterator, N> next_nodes;
    YType operator()(XType X) {
      return (*(*next_nodes[rule(X)]))(X);
    }
    virtual ~RuleNode(void) {}
  };

  struct DataNode : public Node {
    YType y;
    YType operator()(XType X) {
      return y;
    }
    virtual ~DataNode(void) {}
  };

public:
  AbstractTree(void) {}

  virtual ~AbstractTree(void) {
    for (auto node : nodes) {
      delete node;
    }
  }
 
  void fit(std::vector<ExampleType>& data) {
    constructNode(data.begin(), data.end());
  }

  std::vector<YType> predict(std::vector<XType> const& X) {
    std::vector<YType> res;
    res.reserve(X.size());
    for (auto el : X) {
      res.push_back((*(*nodes.begin()))(el));
    }
    return res;
  }

protected:
  std::list<Node*> nodes;


  virtual std::pair<RuleType, bool> generateRule(
    typename std::vector<ExampleType>::const_iterator begin, 
    typename std::vector<ExampleType>::const_iterator end
  ) = 0;

  virtual YType    generateY(
    typename std::vector<ExampleType>::const_iterator begin, 
    typename std::vector<ExampleType>::const_iterator end
  ) = 0;

private:
  std::array<typename std::vector<ExampleType>::iterator, N + 1> segmentData(
    typename std::vector<ExampleType>::iterator begin, 
    typename std::vector<ExampleType>::iterator end, 
    RuleType rule
  ) {
    std::array<typename std::vector<ExampleType>::iterator, N + 1> res;
    res[0] = begin;
    for (size_t i = 0; i < N - 1; ++i) {
      res[i + 1] = std::partition(res[i], end, [i, &rule](auto x) { return rule(x.first) == i; });
    }
    res[N] = end;
    return res;
  }

  typename std::list<ExampleType>::iterator constructNode(
    typename std::vector<ExampleType>::iterator begin, 
    typename std::vector<ExampleType>::iterator end
  ) {
    auto& [rule, type] = generateRule(begin, end);
    Node* node = (type) ? new DataNode(generateY(begin, end)) : new RuleNode(rule);
    nodes.push_back(node);
    auto res = std::prev(nodes.end());
    if (!type) {
      auto& points = segmentData(begin, end, rule);
      for (size_t i = 0; i < N; ++i) {
        node->next_nodes[i] = constructNode(points[i], points[i + 1]);
      }
    }
    return res;
  }

};

} // kaspar
