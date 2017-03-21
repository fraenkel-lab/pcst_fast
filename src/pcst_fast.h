#ifndef __PCST_FAST_H__
#define __PCST_FAST_H__

#include <string>
#include <utility>
#include <vector>

#include "pairing_heap.h"
#include "priority_queue.h"

namespace cluster_approx {

class PCSTFast {
 public:
  enum PruningMethod {
    kNoPruning = 0,
    kSimplePruning,
    kGWPruning,
    kStrongPruning,
    kUnknownPruning,
  };

  struct Statistics {
    long long total_num_edge_events;
    long long num_deleted_edge_events;
    long long num_merged_edge_events;
    long long total_num_merge_events;
    long long num_active_active_merge_events;
    long long num_active_inactive_merge_events;
    long long total_num_edge_growth_events;
    long long num_active_active_edge_growth_events;
    long long num_active_inactive_edge_growth_events;
    long long num_cluster_events;

    Statistics();
  };

  const static int kNoRoot = -1;

  static PruningMethod parse_pruning_method(const std::string& input);


  PCSTFast(const std::vector<std::pair<int, int> >& edges_,
           const std::vector<double>& prizes_,
           const std::vector<double>& costs_,
           int root_,
           int target_num_active_clusters_,
           PruningMethod pruning_,
           int verbosity_level_,
           void (*output_function_)(const char*));
  
  ~PCSTFast();

  bool run(std::vector<int>* result_nodes,
           std::vector<int>* result_edges);

  void get_statistics(Statistics* s);


 private:
  typedef PairingHeap<double, int> PairingHeapType;
  typedef PriorityQueue<double, int> PriorityQueueType;

  struct EdgeInfo {
    int inactive_merge_event;
  };

  struct EdgePart {
    double next_event_val;
    bool deleted;
    PairingHeapType::ItemHandle heap_node;
  };

  struct InactiveMergeEvent {
    int active_cluster_index;
    int inactive_cluster_index;
    int active_cluster_node;
    int inactive_cluster_node;
  };

  struct Cluster {
    PairingHeapType edge_parts;
    bool active;
    double active_start_time;
    double active_end_time;
    int merged_into;
    double prize_sum;
    double subcluster_moat_sum;
    double moat;
    bool contains_root;
    int skip_up;
    double skip_up_sum;
    int merged_along;
    int child_cluster_1;
    int child_cluster_2;
    bool necessary;

    Cluster(std::vector<PairingHeapType::ItemHandle>* heap_buffer)
        : edge_parts(heap_buffer) {}
  };


  const std::vector<std::pair<int, int> >& edges;
  const std::vector<double>& prizes;
  const std::vector<double>& costs;
  int root;
  int target_num_active_clusters;
  PruningMethod pruning;
  int verbosity_level;
  void (*output_function)(const char*);
  Statistics stats;

  std::vector<PairingHeapType::ItemHandle> pairing_heap_buffer;
  std::vector<EdgePart> edge_parts;
  std::vector<EdgeInfo> edge_info;
  std::vector<Cluster> clusters;
  std::vector<InactiveMergeEvent> inactive_merge_events;
  PriorityQueueType clusters_deactivation;
  PriorityQueueType clusters_next_edge_event;
  double current_time;
  double eps;
  std::vector<bool> node_good;  // marks whether a node survives simple pruning
  std::vector<bool> node_deleted;
  std::vector<int> phase2_result;

  std::vector<std::pair<int, double> > path_compression_visited;
  std::vector<int> cluster_queue;
  std::vector<std::vector<std::pair<int, double> > > phase3_neighbors;

  // for strong pruning
  std::vector<int> final_component_label;
  std::vector<std::vector<int> > final_components;
  int root_component_index;
  std::vector<std::pair<int, double> > strong_pruning_parent;
  std::vector<double> strong_pruning_payoff;
  std::vector<std::pair<bool, int> > stack;
  std::vector<int> stack2;
 

  const static int kOutputBufferSize = 10000;
  char output_buffer[kOutputBufferSize];
  
  void get_next_edge_event(double* next_time,
                           int* next_cluster_index,
                           int* next_edge_part_index);

  void remove_next_edge_event(int next_cluster_index);
  
  void get_next_cluster_event(double* next_time, int* next_cluster_index);

  void remove_next_cluster_event();
  
  void get_sum_on_edge_part(int edge_part_index,
                            double* total_sum,
                            double* finished_moat_sum,
                            int* cur_cluser_index);

  void mark_nodes_as_good(int start_cluster_index);

  void mark_clusters_as_necessary(int start_cluster_index);

  void mark_nodes_as_deleted(int start_node_index, int parent_node_index);

  void label_final_component(int start_node_index, int new_component_index);

  void strong_pruning_from(int start_node_index, bool mark_as_deleted);

  int find_best_component_root(int component_index);

  void build_phase1_node_set(const std::vector<int>& edge_set,
                             std::vector<int>* node_set);

  void build_phase3_node_set(std::vector<int>* node_set);

  void build_phase2_node_set(std::vector<int>* node_set);


  int get_other_edge_part_index(int edge_part_index) {
    if (edge_part_index % 2 == 0) {
      return edge_part_index + 1;
    } else {
      return edge_part_index - 1;
    }
  }
};

}  // namespace cluster_approx


#endif
