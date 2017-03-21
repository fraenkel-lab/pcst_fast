pcsf_fast
==============

A library for solving the prize-collecting Steiner forest (PCSF) problem on graphs.

Installation
------------

The core library has no dependencies besides a basic build system for C++11.
Both g++ and clang are currently supported.
The Python wrapper requires a functioning Python build system.

Compile the python wrapper with the supplied makefile:

    make pcst_fast_py

You can then import the package via `import pcst_fast`.

Usage
-----

The `pcst_fast` package contains the following function:

    vertices, edges = pcst_fast(edges, prizes, costs, root, num_clusters, pruning, verbosity_level)
    
The parameters are:
* `edges`: a 2D int64 array. Each row (of length 2) specifies an undirected edge in the input graph. The nodes are labeled 0 to n-1, where n is the number of nodes.
* `prizes`: the node prizes as a 1D float64 array.
* `costs`: the edge costs as a 1D float64 array.
* `root`: the root note for rooted PCST. For the unrooted variant, this parameter should be -1.
* `num_clusters`: the number of connected components in the output.
* `pruning`: a string value indicating the pruning method. Possible values are `'none'`, `'simple'`, `'gw'`, and `'strong'` (all literals are case-insensitive). `'none'` and `'simple'` return intermediate stages of the algorithm and do not have approximation guarantees. They are only intended for development. The standard GW pruning method is `'gw'`, which is also the default. `'strong'` uses "strong pruning", which was introduced in [\[JMP00\]](http://dl.acm.org/citation.cfm?id=338637). It has the same theoretical guarantees as GW pruning but better empirical performance in some cases. For the PCSF problem, the output of strong pruning is at least as good as the output of GW pruning.
* `verbosity_level`: an integer indicating how much debug output the function should produce.

The output variables are:
* `vertices`: the vertices in the solution as a 1D int64 array.
* `edges`: the edges in the output as a 1D int64 array. The list contains indices into the list of edges passed into the function.
