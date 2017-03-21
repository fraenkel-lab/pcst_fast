#include <cstdint>
#include <string>
#include <utility>
#include <vector>

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include "pcst_fast.h"

namespace py = pybind11;

using cluster_approx::PCSTFast;

void output_function(const char* output) {
  py::print(output, py::arg("flush")=true);
}

std::pair<py::array_t<int>, py::array_t<int>> pcst_fast(
    py::array_t<int64_t, py::array::c_style> edges,
    py::array_t<double, py::array::c_style> prizes,
    py::array_t<double, py::array::c_style> costs,
    int root,
    int num_clusters,
    const std::string& pruning,
    int verbosity_level) {
  int target_num_active_clusters = num_clusters;
  if (root >= 0 && num_clusters != 1) {
    throw std::invalid_argument("In the rooted case, only one output cluster "
                                "is supported.");
  }
  if (root >= 0) {
    target_num_active_clusters = 0;
  }
  py::buffer_info edges_info = edges.request();
  if (edges_info.ndim != 2) {
    throw std::invalid_argument("Edges must be a two-dimensional array.");
  }
  if (edges_info.shape[1] != 2) {
    throw std::invalid_argument("The edges array must have two columns.");
  }
  if (edges_info.itemsize != sizeof(int64_t)) {
    throw std::invalid_argument("The edges itemsize does not match "
                                "sizeof(int).");
  }
  int num_edges = edges_info.shape[0];
  int64_t* edges_ptr = static_cast<int64_t*>(edges_info.ptr);
  std::vector<std::pair<int, int> > tmp_edges(num_edges);
  for (int ii = 0; ii < num_edges; ++ii) {
    tmp_edges[ii].first = edges_ptr[2 * ii];
    tmp_edges[ii].second = edges_ptr[2 * ii + 1];
  }

  py::buffer_info prizes_info = prizes.request();
  if (prizes_info.ndim != 1) {
    throw std::invalid_argument("Prizes must be a one-dimensional array.");
  }
  if (prizes_info.itemsize != sizeof(double)) {
    throw std::invalid_argument("The prizes itemsize does not match "
                                "sizeof(double).");
  }
  int num_nodes = prizes_info.shape[0];
  double* prizes_ptr = static_cast<double*>(prizes_info.ptr);
  std::vector<double> tmp_prizes(num_nodes);
  for (int ii = 0; ii < num_nodes; ++ii) {
    tmp_prizes[ii] = prizes_ptr[ii];
  }
  
  py::buffer_info costs_info = costs.request();
  if (costs_info.ndim != 1) {
    throw std::invalid_argument("Costs must be a one-dimensional array.");
  }
  if (costs_info.itemsize != sizeof(double)) {
    throw std::invalid_argument("The costs itemsize does not match "
                                "sizeof(double).");
  }
  if (static_cast<int>(costs_info.shape[0]) != num_edges) {
    throw std::invalid_argument("The size of the costs array does not match "
                                "the number of rows in the edges array.");
  }
  double* costs_ptr = static_cast<double*>(costs_info.ptr);
  std::vector<double> tmp_costs(num_edges);
  for (int ii = 0; ii < num_edges; ++ii) {
    tmp_costs[ii] = costs_ptr[ii];
  }
 
  // TODO: change the code so that fewer copying happens.
  // First, update pcst_fast on the C++ side to use a given memory location
  // as input (instead of a std::vector reference).
  // Similarly, pass in raw C pointers for output.
  
  PCSTFast::PruningMethod pruning_method =
      PCSTFast::parse_pruning_method(pruning);
  PCSTFast algo(tmp_edges, tmp_prizes, tmp_costs, root,
                target_num_active_clusters, pruning_method, verbosity_level,
                output_function);
  std::vector<int> result_nodes;
  std::vector<int> result_edges;
  algo.run(&result_nodes, &result_edges);

  py::array_t<int64_t> result_nodes_array(result_nodes.size());
  py::buffer_info result_nodes_info = result_nodes_array.request();
  int64_t* result_nodes_ptr = static_cast<int64_t*>(result_nodes_info.ptr);
  for (int ii = 0; ii < static_cast<int>(result_nodes.size()); ++ii) {
    result_nodes_ptr[ii] = result_nodes[ii];
  }
  
  py::array_t<int64_t> result_edges_array(result_edges.size());
  py::buffer_info result_edges_info = result_edges_array.request();
  int64_t* result_edges_ptr = static_cast<int64_t*>(result_edges_info.ptr);
  for (int ii = 0; ii < static_cast<int>(result_edges.size()); ++ii) {
    result_edges_ptr[ii] = result_edges[ii];
  }
  
  return std::make_pair(result_nodes_array, result_edges_array);
}

PYBIND11_PLUGIN(pcst_fast) {
  py::module m("pcst_fast", "A fast algorithm for the PCSF problem.");

  m.def("pcst_fast", &pcst_fast, "Runs the pcst_fast algorithm.");

  return m.ptr();
}
