# A junction tree.
struct JunctionTree
    order::Order
    tree::Tree
    seperatorlist::Vector{Vector{Int}}
    residuallist::Vector{Vector{Int}}
    subtreelist::Vector{Int}
    width::Int 
end


# Construct a tree decomposition.
function JunctionTree(graph::AbstractGraph, stree::EliminationTree)
    graph = makeeliminationgraph(graph, stree) 

    m = nv(graph)
    n = length(stree)
    seperatorlist = Vector{Vector{Int}}(undef, n) 
    residuallist = Vector{Vector{Int}}(undef, n)   
 
    seperatorlist[n] = []
    residuallist[n] = stree.firstsupernodelist[n]:m

    for i in 1:n - 1
        v₁ = stree.firstsupernodelist[i]
        v₂ = stree.lastsupernodelist[i]
        v₃ = stree.firstancestorlist[i]
        bag = neighbors(graph, v₁)

        residuallist[i] = v₁:v₂
        seperatorlist[i] = bag[searchsortedfirst(bag, v₃):end]
    end

    order = stree.order
    tree = stree.tree
    subtreelist = stree.subtreelist
    width = stree.width

    JunctionTree(order, tree, seperatorlist, residuallist, subtreelist, width)
end


# Reorient a juncton tree toward a given root.
function JunctionTree(jtree::JunctionTree, root::Integer)
    i = root
    j = parentindex(jtree, i)
    parentlist = copy(jtree.tree.parentlist)
    childrenlist = deepcopy(jtree.tree.childrenlist)

    while !isnothing(j)
        parentlist[j] = i
        push!(childrenlist[i], j)
        deletesorted!(childrenlist[j], i)
        i = j
        j = parentindex(jtree, i)
    end

    tree = Tree(root, parentlist, childrenlist)
    postorder, tree = makepostorder(tree)
    m = length(jtree.subtreelist)
    n = length(jtree)
    seperatorlist = Vector{Vector{Int}}(undef, n)
    residuallist = Vector{Vector{Int}}(undef, n)
    subtreelist = Vector{Int}(undef, m)

    seperatorlist[n] = []
    residuallist[n] = [jtree.residuallist[postorder[n]]; jtree.seperatorlist[postorder[n]]]
    subtreelist[residuallist[n]] .= n    

    for i in 1:n - 1
        j = postorder[i]

        if isdescendant(jtree, root, j)
            i′ = parentindex(tree, i)
            j′ = postorder[i′]
            seperatorlist[i] = jtree.seperatorlist[j′]
            residuallist[i] = [jtree.residuallist[j]; jtree.seperatorlist[j]]
            deletesorted!(residuallist[i], seperatorlist[i])
        else
            seperatorlist[i] = jtree.seperatorlist[j]
            residuallist[i] = jtree.residuallist[j]
        end

        subtreelist[residuallist[i]] .= i
    end 

    order = jtree.order
    width = jtree.width
    JunctionTree(order, tree, seperatorlist, residuallist, subtreelist, width)
end


# Construct a tree decomposition, first computing a supernodal elimination tree.
function JunctionTree(
    graph::AbstractGraph,
    order::Order,
    supernode::Supernode=DEFAULT_SUPERNODE)

    stree = EliminationTree(graph, order, supernode)
    JunctionTree(graph, stree)
end


# Construct a tree decomposition, first computing an elimination order and a supernodal
# elimination tree.
function JunctionTree(
    graph::AbstractGraph,
    algorithm::EliminationAlgorithm=DEFAULT_ELIMINATION_ALGORITHM,
    supernode::Supernode=DEFAULT_SUPERNODE)

    stree = EliminationTree(graph, algorithm, supernode)
    JunctionTree(graph, stree)
end


# Draw a junction tree.
function GraphPlot.gplot(
    jtree::JunctionTree;
    edgelabel::Bool=false,
    nodelabel::Bool=false,
    nodefillc=RGB(0.1,0.1,0.1),
    rootfillc=RGB(0.882911,0.359638,0.360092),
    kwargs...)

    el = edgelabel
    nl = nodelabel
    nc = parse(RGB{Float64}, nodefillc)
    rc = parse(RGB{Float64}, rootfillc)

    n = length(jtree)
    graph = Graph(n)
    edgelabel = Vector{Vector{Int}}(undef, n - 1)
    nodelabel = Vector{Vector{Int}}(undef, n)
    nodefillc = Vector{RGB{Float64}}(undef, n)

    nodelabel[n] = getresidual(jtree, n)
    nodefillc[n] = rc 

    for i in n - 1:-1:1
        add_edge!(graph, i, parentindex(jtree, i))
        edgelabel[i] = getseperator(jtree, i)
        nodelabel[i] = getresidual(jtree, i)
        nodefillc[i] = nc
    end

    nodelabel = nl ? nodelabel : nothing
    edgelabel = el ? edgelabel : []
    nodelabelc = "white"

    gplot(graph; edgelabel, nodelabel, nodefillc, nodelabelc, kwargs...)
end


function Base.length(jtree::JunctionTree)
    length(jtree.tree)
end


function getwidth(jtree::JunctionTree)
    jtree.width
end


# Get the seperator at node i.
function getseperator(jtree::JunctionTree, i::Integer)
    jtree.order[jtree.seperatorlist[i]]
end


# Get the residual at node i.
function getresidual(jtree::JunctionTree, i::Integer)
    jtree.order[jtree.residuallist[i]]
end


# Get the highest node containing the vertex v.
function getsubtree(jtree::JunctionTree, v::Integer)
    jtree.subtreelist[jtree.order.index[v]]
end

# Get the highest node containing the vertices vs.
function getsubtree(jtree::JunctionTree, vs::AbstractVector)
    nodes = map(vs) do v
        getsubtree(jtree, v)
    end

    minimum(nodes)
end


# Get the level of node i.
function getlevel(jtree::JunctionTree, i::Integer)
    getlevel(jtree.tree, i)
end


# Evaluate whether node i₁ is a descendant of node i₂.
function AbstractTrees.isdescendant(jtree::JunctionTree, i₁::Integer, i₂::Integer)
    isdescendant(jtree.tree, i₁, i₂)
end


# Construct an elimination graph.
function makeeliminationgraph(graph::AbstractGraph, stree::EliminationTree)
    n = length(stree)
    graph = DiGraph(graph, stree.order)

    for i in 1:n - 1
        u₁ = stree.firstsupernodelist[i]
        u₂ = stree.lastsupernodelist[i]

        for u in u₁:u₂ - 1
            v = u + 1

            for w in neighbors(graph, u)
                if v != w
                    add_edge!(graph, v, w)
                end
            end
        end

        u = u₂
        v = stree.firstsupernodelist[parentindex(stree, i)]

        for w in neighbors(graph, u)
            if v != w
                add_edge!(graph, v, w)
            end
        end
    end

    graph
end


##########################
# Indexed Tree Interface #
##########################


function AbstractTrees.rootindex(jtree::JunctionTree)
    rootindex(jtree.tree)
end


function AbstractTrees.parentindex(jtree::JunctionTree, i::Integer)
    parentindex(jtree.tree, i)
end


function AbstractTrees.childindices(jtree::JunctionTree, i::Integer)
    childindices(jtree.tree, i)
end


function AbstractTrees.NodeType(::Type{IndexNode{JunctionTree, Int}})
    HasNodeType()
end


function AbstractTrees.nodetype(::Type{IndexNode{JunctionTree, Int}})
    IndexNode{JunctionTree, Int}
end
