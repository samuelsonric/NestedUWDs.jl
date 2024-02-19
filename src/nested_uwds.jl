struct NestedUWD{T <: UndirectedWiringDiagram}
    diagram::T
    jtree::JunctionTree
    assignments::Vector{Int}
end


function NestedUWD(
    diagram::UndirectedWiringDiagram,
    hypergraph::Hypergraph,
    jtree::JunctionTree)

    m = Graphs.ne(hypergraph) - 1
    query = getsubtree(jtree, Graphs.vertices(hypergraph, 1 + m))
    jtree = JunctionTree(jtree, query)

    assignments = map(1:m) do e
        getsubtree(jtree, Graphs.vertices(hypergraph, e))
    end

    NestedUWD(diagram, jtree, assignments)
end 


function NestedUWD(
    diagram::UndirectedWiringDiagram,
    order::Order,
    supernode::Supernode=DEFAULT_SUPERNODE)

    hypergraph = Hypergraph(diagram)
    graph = Graphs.Graph(hypergraph)
    jtree = JunctionTree(graph, order, supernode)
    NestedUWD(diagram, hypergraph, jtree)
end


function NestedUWD(
    diagram::UndirectedWiringDiagram,
    algorithm::EliminationAlgorithm=DEFAULT_ELIMINATION_ALGORITHM,
    supernode::Supernode=DEFAULT_SUPERNODE)

    hypergraph = Hypergraph(diagram)
    graph = Graphs.Graph(hypergraph)
    jtree = JunctionTree(graph, algorithm, supernode)
    NestedUWD(diagram, hypergraph, jtree)
end


# Given an undirected wiring diagram D, construct a hypergraph whose edges are the boxes of
# D, and whose vertices are the junctions of D.
function JunctionTrees.Hypergraph(diagram::UndirectedWiringDiagram)
    m = 1 + nboxes(diagram)
    n = njunctions(diagram)
    e2v = Vector{Vector{Int}}(undef, m)
    v2e = Vector{Vector{Int}}(undef, n)
    v2e .= [[]]
    e2v .= [[]]

    for p in ports(diagram)
        b = box(diagram, p)
        j = junction(diagram, p)
        push!(e2v[b], j)
        push!(v2e[j], b)
    end

    for p in ports(diagram; outer=true)
        j = junction(diagram, p; outer=true)
        push!(e2v[m], j)
        push!(v2e[j], m)
    end

    Hypergraph(e2v, v2e)
end


function makeschedule(nuwd::NestedUWD{<:Union{
    UntypedRelationDiagram{T},
    UntypedHypergraphDiagram{T},
    TypedRelationDiagram{<:Any, T},
    HypergraphDiagram{<:Any, T}}}) where T

    m = length(nuwd.assignments)
    n = length(nuwd.jtree)

    parents = map(1:n - 1) do i
        parentindex(nuwd.jtree, i)
    end

    schedule = WiringDiagramACSet{T, Nothing, Union{Nothing, AbstractBox}, DataType}()

    add_parts!(schedule, :Box, n)
    add_parts!(schedule, :Wire, n - 1)
    add_parts!(schedule, :InPort, m + n - 1)
    add_parts!(schedule, :InWire, m)
    add_parts!(schedule, :OutPort, n)
    add_parts!(schedule, :OutWire, 1)
    add_parts!(schedule, :OuterInPort, m)
    add_parts!(schedule, :OuterOutPort, 1)

    schedule[:, :src] = 1:n - 1
    schedule[:, :tgt] = m + 1:m + n - 1
    schedule[:, :in_src] = 1:m
    schedule[:, :in_tgt] = 1:m
    schedule[:, :out_src] = n:n
    schedule[:, :out_tgt] = 1:1
    schedule[:, :in_port_box] = [nuwd.assignments; parents]
    schedule[:, :out_port_box] = 1:n

    schedule[:, :box_type] = Box{Nothing}
    schedule[:, :outer_in_port_type] = nuwd.diagram[:, :name]

    Theory = ThSymmetricMonoidalCategory.Meta.T
    WiringDiagram{Theory, T, Nothing, Nothing}(schedule, nothing)
end
