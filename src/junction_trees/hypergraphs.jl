# A hypergraph.
struct Hypergraph
    e2v::Vector{Vector{Int}}
    v2e::Vector{Vector{Int}}
end


# Get the edges incident to a vertex v.
function Graphs.edges(hypergraph::Hypergraph, v::Integer)
    hypergraph.v2e[v]
end


# Get the vertices incident to an edge e.
function Graphs.vertices(hypergraph::Hypergraph, e::Integer)
    hypergraph.e2v[e]
end


# Get the edges of a hypergraph.
function Graphs.edges(hypergraph::Hypergraph)
    1:ne(hypergraph)
end


# Get the vertices of a hypergraph.
function Graphs.vertices(hypergraph::Hypergraph)
    1:nv(hypergraph)
end


# Get the number of edges incident to a vertex v.
function Graphs.ne(hypergraph::Hypergraph, v::Integer)
    length(edges(hypergraph, v))
end


# Get the number of vertices incident to an edge e.
function Graphs.nv(hypergraph::Hypergraph, e::Integer)
    length(vertices(hypergraph, e))
end


# Get the number of edges in a hypergraph.
function Graphs.ne(hypergraph::Hypergraph)
    length(hypergraph.e2v)
end


# Get the number of vertices in a hypergraph.
function Graphs.nv(hypergraph::Hypergraph)
    length(hypergraph.v2e)
end


# Construct the primal graph of a hypergraph.
function Graphs.Graph(hypergraph::Hypergraph)
    graph = Graph(nv(hypergraph))

    for e in edges(hypergraph)
        n = nv(hypergraph, e)
        V = vertices(hypergraph, e)
        
        for i in 1:n, j in i + 1:n
            add_edge!(graph, V[i], V[j])
        end
    end

    graph
end
