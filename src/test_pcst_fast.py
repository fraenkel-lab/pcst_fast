import numpy as np
import numpy.testing as npt
import sys

sys.path.append('..')

import pcst_fast

def test_simple1():
  edges = [[0, 1], [1, 2], [2, 3]]
  prizes = [1.0, 1.0, 1.0, 1.0]
  costs = [0.8, 1.8, 2.8]
  result_nodes, result_edges = pcst_fast.pcst_fast(edges, prizes, costs, -1, 1, 'strong', 0)
  npt.assert_array_equal(result_nodes, [0, 1])
  npt.assert_array_equal(result_edges, [0])

def test_simple2():
  edges = [[0, 1], [1, 2]]
  prizes = [0.0, 5.0, 6.0]
  costs = [3.0, 4.0]
  result_nodes, result_edges = pcst_fast.pcst_fast(edges, prizes, costs, 0, 1, 'none', 0)
  npt.assert_array_equal(sorted(result_nodes), [0, 1, 2])
  npt.assert_array_equal(sorted(result_edges), [0, 1])
