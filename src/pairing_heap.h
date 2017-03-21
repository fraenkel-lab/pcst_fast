#ifndef __PAIRING_HEAP_H__
#define __PAIRING_HEAP_H__

// TODO: remove printfs
//#include <cstdio>
#include <cstdlib>
#include <vector>

namespace cluster_approx {

template <typename ValueType, typename PayloadType>
class PairingHeap {
 private:
  struct Node {
    Node* sibling;
    Node* child;
    Node* left_up;
    ValueType value;
    ValueType child_offset;
    PayloadType payload;
  };
 
 public:
  typedef Node* ItemHandle;

  PairingHeap(std::vector<ItemHandle>* shared_buffer) : root(NULL) {
    buffer = shared_buffer;
  }

  void release_memory() {
    // Delete heap nodes
    buffer->resize(0);
    if (root != NULL) {
      buffer->push_back(root);
    }
    size_t curi = 0;
    while (curi < buffer->size()) {
      Node* cur_node = (*buffer)[curi];
      if (cur_node->child != NULL) {
        buffer->push_back(cur_node->child);
      }
      if (cur_node->sibling != NULL) {
        buffer->push_back(cur_node->sibling);
      }
      curi += 1;
    }
    for (size_t ii = 0; ii < buffer->size(); ++ii) {
      delete (*buffer)[ii];
    }
  }

  bool is_empty() {
    return root == NULL;
  }

  bool get_min(ValueType* value, PayloadType* payload) {
    if (root != NULL) {
      //printf("In get_min, current root: %x\n", root);
      //fflush(stdout);
      *value = root->value;
      *payload = root->payload;
      return true;
    } else {
      return false;
    }
  }

  ItemHandle insert(ValueType value, PayloadType payload) {
    Node* new_node = new Node();
    new_node->sibling = NULL;
    new_node->child = NULL;
    new_node->left_up = NULL;
    new_node->value = value;
    new_node->payload = payload;
    new_node->child_offset = 0;
    root = link(root, new_node);
    return new_node;
  }

  void add_to_heap(ValueType value) {
    if (root != NULL) {
      root->value += value;
      root->child_offset += value;
    }
  }

  void decrease_key(ItemHandle node, ValueType from_value, ValueType to_value) {
    ValueType additional_offset = from_value - node->value;
    node->child_offset += additional_offset;
    node->value = to_value;
    if (node->left_up != NULL) {
      if (node->left_up->child == node) {
        node->left_up->child = node->sibling;
      } else {
        node->left_up->sibling = node->sibling;
      }
      if (node->sibling != NULL) {
        node->sibling->left_up = node->left_up;
      }
      node->left_up = NULL;
      node->sibling = NULL;
      root = link(root, node);
    }
  }

  bool delete_min(ValueType* value, PayloadType* payload) {
    if (root == NULL) {
      return false;
    }
    //printf("In delete_min, root is %x (payload %d)\n", root, root->payload);
    //fflush(stdout);
    Node* result = root;
    buffer->resize(0);
    Node* cur_child = root->child;
    Node* next_child = NULL;
    while (cur_child != NULL) {
      buffer->push_back(cur_child);
      //printf("In delete_min, added child %x to buffer\n", cur_child);
      //fflush(stdout);
      next_child = cur_child->sibling;
      cur_child->left_up = NULL;
      cur_child->sibling = NULL;
      cur_child->value += result->child_offset;
      cur_child->child_offset += result->child_offset;
      cur_child = next_child;
    }

    //printf("In delete_min, root hat %lu children\n", buffer->size());
    //fflush(stdout);

    size_t merged_children = 0;
    while (merged_children + 2 <= buffer->size()) {
      (*buffer)[merged_children / 2] = link(
          (*buffer)[merged_children], (*buffer)[merged_children + 1]);
      merged_children += 2;
    }
    if (merged_children != buffer->size()) {
      (*buffer)[merged_children / 2] = (*buffer)[merged_children];
      buffer->resize(merged_children / 2 + 1);
    } else {
      buffer->resize(merged_children / 2);
    }

    if (buffer->size() > 0) {
      root = (*buffer)[buffer->size() - 1];
      for (int ii = buffer->size() - 2; ii >= 0; --ii) {
        root = link(root, (*buffer)[ii]);
      }
    } else {
      root = NULL;
    }

    *value = result->value;
    *payload = result->payload;
    //printf("In delete_min, deleting %x\n", result);
    //printf("In delete_min, new root: %x\n", root);
    //fflush(stdout);
    delete result;
    return true;
  }

  static PairingHeap meld(PairingHeap* heap1, PairingHeap* heap2) {
    PairingHeap result(heap1->buffer);
    result.root = link(heap1->root, heap2->root);
    heap1->root = NULL;
    heap2->root = NULL;
    return result;
  }

 private:
  Node* root;
  std::vector<ItemHandle>* buffer;

  static Node* link(Node* node1, Node* node2) {
    if (node1 == NULL) {
      return node2;
    }
    if (node2 == NULL) {
      return node1;
    }
    Node* smaller_node = node2;
    Node* larger_node = node1;
    if (node1->value < node2->value) {
      smaller_node = node1;
      larger_node = node2;
    }
    //printf("Linking %x (smaller node) and %x (larger node)\n",
    //    smaller_node, larger_node);
    //fflush(stdout);
    larger_node->sibling = smaller_node->child;
    if (larger_node->sibling != NULL) {
      larger_node->sibling->left_up = larger_node;
    }
    larger_node->left_up = smaller_node;
    smaller_node->child = larger_node;
    larger_node->value -= smaller_node->child_offset;
    larger_node->child_offset -= smaller_node->child_offset;
    return smaller_node;
  }
};

}  // namespace cluster_approx

#endif
