#include "pcst_fast.h"

#include <algorithm>
#include <utility>
#include <vector>

#include "gtest/gtest.h"

#include "test_helpers.h"

using cluster_approx::PCSTFast;
using std::make_pair;
using std::pair;
using std::vector;

const int kVerbosityLevel = 0;


void RunAlgo(const vector<pair<int, int> >& edges,
             const vector<double>& prizes,
             const vector<double>& costs,
             int root,
             int target_num_active_clusters,
             PCSTFast::PruningMethod pruning,
             const vector<int>& expected_node_result,
             const vector<int>& expected_edge_result) {

  vector<int> node_result;
  vector<int> edge_result;
  PCSTFast algo(edges, prizes, costs, root, target_num_active_clusters,
                pruning, kVerbosityLevel, WriteToStderr);
  ASSERT_TRUE(algo.run(&node_result, &edge_result));

  std::sort(node_result.begin(), node_result.end());
  std::sort(edge_result.begin(), edge_result.end());
  vector<int> sorted_expected_node_result(expected_node_result);
  std::sort(sorted_expected_node_result.begin(),
            sorted_expected_node_result.end());
  vector<int> sorted_expected_edge_result(expected_edge_result);
  std::sort(sorted_expected_edge_result.begin(),
            sorted_expected_edge_result.end());

  CheckResult(sorted_expected_node_result, node_result);
  CheckResult(sorted_expected_edge_result, edge_result);
}


template <size_t N1, size_t N2, size_t N3, size_t N4>
void RunAlgo(const vector<pair<int, int> >& edges,
             const double (&prizes)[N1],
             const double (&costs)[N2],
             int root,
             int target_num_active_clusters,
             PCSTFast::PruningMethod pruning,
             const int (&expected_node_result)[N3],
             const int (&expected_edge_result)[N4]) {
  vector<double> prizes_(begin(prizes), end(prizes));
  vector<double> costs_(begin(costs), end(costs));
  vector<int> expected_node_result_(begin(expected_node_result),
                                    end(expected_node_result));
  vector<int> expected_edge_result_(begin(expected_edge_result),
                                    end(expected_edge_result));
  RunAlgo(edges, prizes_, costs_, root, target_num_active_clusters, pruning,
          expected_node_result_, expected_edge_result_);
}


template <size_t N1, size_t N2, size_t N3>
void RunAlgo(const vector<pair<int, int> >& edges,
             const double (&prizes)[N1],
             const double (&costs)[N2],
             int root,
             int target_num_active_clusters,
             PCSTFast::PruningMethod pruning,
             const int (&expected_node_result)[N3]) {
  vector<double> prizes_(begin(prizes), end(prizes));
  vector<double> costs_(begin(costs), end(costs));
  vector<int> expected_node_result_(begin(expected_node_result),
                                    end(expected_node_result));
  vector<int> expected_edge_result_(0);
  RunAlgo(edges, prizes_, costs_, root, target_num_active_clusters, pruning,
          expected_node_result_, expected_edge_result_);
}


TEST(PCSTFastTest, SimpleTestRootedNoPruning) {
  vector<pair<int, int> > edges;
  edges.push_back(make_pair(0, 1));
  edges.push_back(make_pair(1, 2));
  const double prizes[] = {0, 5, 6};
  const double costs[] = {3, 4};
  int root = 0;
  int target_num_active_clusters = 0;
  PCSTFast::PruningMethod pruning = PCSTFast::kNoPruning;

  const int node_result[] = {0, 1, 2};
  const int edge_result[] = {0, 1};

  RunAlgo(edges, prizes, costs, root, target_num_active_clusters, pruning,
          node_result, edge_result);
}


TEST(PCSTFastTest, SimpleTestUnrootedNoPruning) {
  vector<pair<int, int> > edges;
  edges.push_back(make_pair(0, 1));
  edges.push_back(make_pair(1, 2));
  const double prizes[] = {0, 5, 6};
  const double costs[] = {3, 4};
  int root = -1;
  int target_num_active_clusters = 1;
  PCSTFast::PruningMethod pruning = PCSTFast::kNoPruning;

  const int node_result[] = {1, 2};
  const int edge_result[] = {1};

  RunAlgo(edges, prizes, costs, root, target_num_active_clusters, pruning,
          node_result, edge_result);
}


TEST(PCSTFastTest, SimpleTestUnrootedGWPruning) {
  vector<pair<int, int> > edges;
  edges.push_back(make_pair(0, 1));
  edges.push_back(make_pair(1, 2));
  const double prizes[] = {0, 5, 6};
  const double costs[] = {3, 4};
  int root = -1;
  int target_num_active_clusters = 1;
  PCSTFast::PruningMethod pruning = PCSTFast::kGWPruning;

  const int node_result[] = {1, 2};
  const int edge_result[] = {1};

  RunAlgo(edges, prizes, costs, root, target_num_active_clusters, pruning,
          node_result, edge_result);
}


TEST(PCSTFastTest, SimpleTestUnrootedStrongPruning) {
  vector<pair<int, int> > edges;
  edges.push_back(make_pair(0, 1));
  edges.push_back(make_pair(1, 2));
  const double prizes[] = {0, 5, 6};
  const double costs[] = {3, 4};
  int root = -1;
  int target_num_active_clusters = 1;
  PCSTFast::PruningMethod pruning = PCSTFast::kStrongPruning;

  const int node_result[] = {1, 2};
  const int edge_result[] = {1};

  RunAlgo(edges, prizes, costs, root, target_num_active_clusters, pruning,
          node_result, edge_result);
}


TEST(PCSTFastTest, Simple2TestRootedNoPruning) {
  vector<pair<int, int> > edges;
  edges.push_back(make_pair(0, 1));
  edges.push_back(make_pair(1, 2));
  edges.push_back(make_pair(2, 3));
  const double prizes[] = {10, 0, 1, 10};
  const double costs[] = {10, 4, 3};
  int root = 0;
  int target_num_active_clusters = 0;
  PCSTFast::PruningMethod pruning = PCSTFast::kNoPruning;

  const int node_result[] = {0, 1, 2, 3};
  const int edge_result[] = {1, 2};

  RunAlgo(edges, prizes, costs, root, target_num_active_clusters, pruning,
          node_result, edge_result);
}


TEST(PCSTFastTest, Simple2TestRootedGWPruning) {
  vector<pair<int, int> > edges;
  edges.push_back(make_pair(0, 1));
  edges.push_back(make_pair(1, 2));
  edges.push_back(make_pair(2, 3));
  const double prizes[] = {10, 0, 1, 10};
  const double costs[] = {10, 4, 3};
  int root = 0;
  int target_num_active_clusters = 0;
  PCSTFast::PruningMethod pruning = PCSTFast::kGWPruning;

  const int node_result[] = {0};
  // edge result should be empty

  RunAlgo(edges, prizes, costs, root, target_num_active_clusters, pruning,
          node_result);
}


TEST(PCSTFastTest, Simple3TestRootedNoPruning) {
  vector<pair<int, int> > edges;
  edges.push_back(make_pair(0, 1));
  edges.push_back(make_pair(1, 2));
  edges.push_back(make_pair(2, 3));
  const double prizes[] = {10, 10, 1, 10};
  const double costs[] = {10, 6, 5};
  int root = 0;
  int target_num_active_clusters = 0;
  PCSTFast::PruningMethod pruning = PCSTFast::kNoPruning;

  const int node_result[] = {0, 1, 2, 3};
  const int edge_result[] = {0, 1, 2};

  RunAlgo(edges, prizes, costs, root, target_num_active_clusters, pruning,
          node_result, edge_result);
}


TEST(PCSTFastTest, Simple3TestRootedGWPruning) {
  vector<pair<int, int> > edges;
  edges.push_back(make_pair(0, 1));
  edges.push_back(make_pair(1, 2));
  edges.push_back(make_pair(2, 3));
  const double prizes[] = {10, 10, 1, 10};
  const double costs[] = {10, 6, 5};
  int root = 0;
  int target_num_active_clusters = 0;
  PCSTFast::PruningMethod pruning = PCSTFast::kGWPruning;

  const int node_result[] = {0, 1, 2, 3};
  const int edge_result[] = {0, 1, 2};

  RunAlgo(edges, prizes, costs, root, target_num_active_clusters, pruning,
          node_result, edge_result);
}


TEST(PCSTFastTest, Simple4TestRootedNoPruning) {
  vector<pair<int, int> > edges;
  edges.push_back(make_pair(0, 1));
  edges.push_back(make_pair(1, 2));
  const double prizes[] = {10, 3, 3};
  const double costs[] = {100, 2};
  int root = 0;
  int target_num_active_clusters = 0;
  PCSTFast::PruningMethod pruning = PCSTFast::kNoPruning;

  const int node_result[] = {0, 1, 2};
  const int edge_result[] = {1};

  RunAlgo(edges, prizes, costs, root, target_num_active_clusters, pruning,
          node_result, edge_result);
}


TEST(PCSTFastTest, Simple4TestRootedGWPruning) {
  vector<pair<int, int> > edges;
  edges.push_back(make_pair(0, 1));
  edges.push_back(make_pair(1, 2));
  const double prizes[] = {10, 3, 3};
  const double costs[] = {100, 2};
  int root = 0;
  int target_num_active_clusters = 0;
  PCSTFast::PruningMethod pruning = PCSTFast::kGWPruning;

  const int node_result[] = {0};
  // no edges expected

  RunAlgo(edges, prizes, costs, root, target_num_active_clusters, pruning,
          node_result);
}


TEST(PCSTFastTest, Simple4TestUnRootedGWPruning) {
  vector<pair<int, int> > edges;
  edges.push_back(make_pair(0, 1));
  edges.push_back(make_pair(1, 2));
  const double prizes[] = {10, 3, 3};
  const double costs[] = {100, 2};
  int root = -1;
  int target_num_active_clusters = 2;
  PCSTFast::PruningMethod pruning = PCSTFast::kGWPruning;

  const int node_result[] = {0, 1, 2};
  const int edge_result[] = {1};

  RunAlgo(edges, prizes, costs, root, target_num_active_clusters, pruning,
          node_result, edge_result);
}


TEST(PCSTFastTest, Simple4bTestUnRootedGWPruning) {
  vector<pair<int, int> > edges;
  edges.push_back(make_pair(0, 1));
  edges.push_back(make_pair(1, 2));
  const double prizes[] = {10, 3, 3};
  const double costs[] = {100, 2};
  int root = -1;
  int target_num_active_clusters = 1;
  PCSTFast::PruningMethod pruning = PCSTFast::kGWPruning;

  const int node_result[] = {0};
  // no edges expected

  RunAlgo(edges, prizes, costs, root, target_num_active_clusters, pruning,
          node_result);
}


TEST(PCSTFastTest, Simple5TestUnRootedGWPruning) {
  vector<pair<int, int> > edges;
  edges.push_back(make_pair(0, 1));
  edges.push_back(make_pair(1, 2));
  edges.push_back(make_pair(2, 3));
  const double prizes[] = {10, 0, 6, 6};
  const double costs[] = {100, 2, 5};
  int root = -1;
  int target_num_active_clusters = 2;
  PCSTFast::PruningMethod pruning = PCSTFast::kGWPruning;

  const int node_result[] = {0, 2, 3};
  const int edge_result[] = {2};

  RunAlgo(edges, prizes, costs, root, target_num_active_clusters, pruning,
          node_result, edge_result);
}


TEST(PCSTFastTest, Medium1TestRootedGWPruning) {
  vector<pair<int, int> > edges;
  edges.push_back(make_pair(0, 1));
  edges.push_back(make_pair(1, 2));
  edges.push_back(make_pair(2, 3));
  edges.push_back(make_pair(0, 9));
  edges.push_back(make_pair(0, 2));
  edges.push_back(make_pair(0, 3));
  edges.push_back(make_pair(0, 5));
  edges.push_back(make_pair(1, 9));
  edges.push_back(make_pair(1, 3));
  edges.push_back(make_pair(1, 5));
  edges.push_back(make_pair(1, 7));
  edges.push_back(make_pair(2, 8));
  edges.push_back(make_pair(2, 3));
  edges.push_back(make_pair(3, 4));
  edges.push_back(make_pair(3, 5));
  edges.push_back(make_pair(3, 6));
  edges.push_back(make_pair(3, 7));
  edges.push_back(make_pair(3, 8));
  edges.push_back(make_pair(3, 9));
  edges.push_back(make_pair(4, 5));
  edges.push_back(make_pair(4, 6));
  edges.push_back(make_pair(4, 7));
  edges.push_back(make_pair(5, 8));
  edges.push_back(make_pair(6, 8));
  const double prizes[] = {0.032052554364677466,
                           0.32473378289799926,
                           0.069699345546302638,
                           0,
                           0.74867253235151754,
                           0.19804330340026255,
                           0.85430521133171622,
                           0.83819939651391351,
                           0.71744625276884877,
                           0.016798567754083948};
  const double costs[] = {0.8,
                          0.8,
                          0.8800000000000001,
                          0.8,
                          0.8,
                          0.8800000000000001,
                          0.8,
                          0.8,
                          0.8800000000000001,
                          0.8,
                          0.8,
                          0.8,
                          0.8800000000000001,
                          0.8800000000000001,
                          0.8800000000000001,
                          0.8800000000000001,
                          0.8800000000000001,
                          0.8800000000000001,
                          0.8800000000000001,
                          0.8,
                          0.8,
                          0.8,
                          0.8,
                          0.8};
  int root = 3;
  int target_num_active_clusters = 0;
  PCSTFast::PruningMethod pruning = PCSTFast::kGWPruning;

  const int node_result[] = {3, 4, 6, 7, 8};
  const int edge_result[] = {16, 20, 21, 23};

  RunAlgo(edges, prizes, costs, root, target_num_active_clusters, pruning,
          node_result, edge_result);
}


TEST(PCSTFastTest, Simple6TestUnRootedGWPruning) {
  vector<pair<int, int> > edges;
  edges.push_back(make_pair(0, 1));
  edges.push_back(make_pair(1, 2));
  edges.push_back(make_pair(2, 3));
  edges.push_back(make_pair(3, 4));
  edges.push_back(make_pair(4, 5));
  edges.push_back(make_pair(5, 6));
  edges.push_back(make_pair(6, 7));
  const double prizes[] = {100.0,
                           0.0,
                           0.0,
                           1.0,
                           0.0,
                           0.0,
                           0.0,
                           100.0};
  const double costs[] = {0.9,
                          0.9,
                          0.9,
                          0.9,
                          0.9,
                          0.9,
                          0.9};
  int root = -1;
  int target_num_active_clusters = 1;
  PCSTFast::PruningMethod pruning = PCSTFast::kGWPruning;

  const int node_result[] = {0, 1, 2, 3, 4, 5, 6, 7};
  const int edge_result[] = {0, 1, 2, 3, 4, 5, 6};

  RunAlgo(edges, prizes, costs, root, target_num_active_clusters, pruning,
          node_result, edge_result);
}


TEST(PCSTFastTest, Simple7TestUnrootedStrongPruning) {
  vector<pair<int, int> > edges;
  edges.push_back(make_pair(0, 1));
  edges.push_back(make_pair(0, 2));
  edges.push_back(make_pair(2, 3));
  edges.push_back(make_pair(3, 4));
  const double prizes[] = {0, 2.2, 0, 0, 2.1};
  const double costs[] = {1, 1, 1, 1};
  int root = -1;
  int target_num_active_clusters = 1;
  PCSTFast::PruningMethod pruning = PCSTFast::kStrongPruning;

  const int node_result[] = {1};
  // no edges expected

  RunAlgo(edges, prizes, costs, root, target_num_active_clusters, pruning,
          node_result);
}


TEST(PCSTFastTest, Simple7TestUnrootedGWPruning) {
  vector<pair<int, int> > edges;
  edges.push_back(make_pair(0, 1));
  edges.push_back(make_pair(0, 2));
  edges.push_back(make_pair(2, 3));
  edges.push_back(make_pair(3, 4));
  const double prizes[] = {0, 2.2, 0, 0, 2.1};
  const double costs[] = {1, 1, 1, 1};
  int root = -1;
  int target_num_active_clusters = 1;
  PCSTFast::PruningMethod pruning = PCSTFast::kGWPruning;

  const int node_result[] = {0, 1, 2, 3, 4};
  const int edge_result[] = {0, 1, 2, 3};

  RunAlgo(edges, prizes, costs, root, target_num_active_clusters, pruning,
          node_result, edge_result);
}


TEST(PCSTFastTest, Simple8TestUnrootedStrongPruning) {
  vector<pair<int, int> > edges;
  edges.push_back(make_pair(0, 1));
  edges.push_back(make_pair(1, 2));
  const double prizes[] = {2, 2, 2};
  const double costs[] = {0, 5};
  int root = -1;
  int target_num_active_clusters = 1;
  PCSTFast::PruningMethod pruning = PCSTFast::kStrongPruning;

  const int node_result[] = {0, 1};
  const int edge_result[] = {0};

  RunAlgo(edges, prizes, costs, root, target_num_active_clusters, pruning,
          node_result, edge_result);
}
