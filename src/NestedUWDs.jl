module NestedUWDs

import Graphs

using AbstractTrees
using Catlab.ACSetInterface
using Catlab.DirectedWiringDiagrams
using Catlab.DirectedWiringDiagrams: WiringDiagramACSet
using Catlab.MonoidalUndirectedWiringDiagrams
using Catlab.MonoidalUndirectedWiringDiagrams: UntypedHypergraphDiagram
using Catlab.RelationalPrograms
using Catlab.Theories
using Catlab.UndirectedWiringDiagrams

include("./junction_trees/JunctionTrees.jl")
using .JunctionTrees
using .JunctionTrees: Hypergraph, DEFAULT_ELIMINATION_ALGORITHM, DEFAULT_SUPERNODE

# Elimination Algorithms
export EliminationAlgorithm, AMDJL_AMD, CuthillMcKeeJL_RCM, MetisJL_ND, MCS

# Supernodes
export Supernode, Node, MaximalSupernode, FundamentalSupernode

# Elimination Trees
export EliminationTree
export getwidth, getsupernode, getsubtree, getlevel

# Junction Trees
export JunctionTree
export getseperator, getresidual

# Nested UWDs
export NestedUWD
export makeschedule

include("./nested_uwds.jl")

end
