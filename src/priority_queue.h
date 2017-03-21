#ifndef __PRIORITY_QUEUE_H__
#define __PRIORITY_QUEUE_H__

#include <map>
#include <set>
#include <vector>

namespace cluster_approx {

template <typename ValueType, typename IndexType>
class PriorityQueue {
 public:
  PriorityQueue() {}

  bool is_empty() {
    return sorted_set.empty();
  }

  bool get_min(ValueType* value, IndexType* index) {
    if (sorted_set.empty()) {
      return false;
    }
    *value = sorted_set.begin()->first;
    *index = sorted_set.begin()->second;
    return true;
  }

  bool delete_min(ValueType* value, IndexType* index) {
    if (sorted_set.empty()) {
      return false;
    }
    *value = sorted_set.begin()->first;
    *index = sorted_set.begin()->second;
    sorted_set.erase(sorted_set.begin());
    return true;
  }

  void insert(ValueType value, IndexType index) {
    if (index >= static_cast<int>(index_to_iterator.size())) {
      index_to_iterator.resize(index + 1);
    }
    index_to_iterator[index] =
        sorted_set.insert(std::make_pair(value, index)).first;
  }

  void decrease_key(ValueType new_value, IndexType index) {
    sorted_set.erase(index_to_iterator[index]);
    index_to_iterator[index] =
        sorted_set.insert(std::make_pair(new_value, index)).first;
  }

  void delete_element(IndexType index) {
    sorted_set.erase(index_to_iterator[index]);
  }
 
 private:
  std::set<std::pair<ValueType, IndexType> > sorted_set;
  std::vector<typename std::set<std::pair<ValueType, IndexType> >::iterator>
      index_to_iterator;
};

}  // namespace cluster_approx

#endif
