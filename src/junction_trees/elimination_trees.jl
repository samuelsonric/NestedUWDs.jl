# A supernodal elimination tree.
struct EliminationTree
    order::Order
    tree::Tree
    firstancestorlist::Vector{Int}
    firstsupernodelist::Vector{Int}
    lastsupernodelist::Vector{Int}
    subtreelist::Vector{Int}
    width::Int
end


# Construct a supernodal elimination tree.
#
# The complexity is
# 𝒪(m α(m, n) + n)
# where m = |E|, n = |V|, and α is the inverse Ackermann function.
function EliminationTree(
    graph::AbstractGraph,
    order::Order,
    supernode::Supernode=DEFAULT_SUPERNODE)

    etree = Tree(graph, order)
    _, outdegreelist = getdegrees(graph, order, etree)

    supernodelist, subtreelist, parentlist, firstancestorlist = makestree(
        etree,
        outdegreelist,
        supernode)

    n = nv(graph)
    tree = Tree(subtreelist[n], parentlist)

    m = length(tree)
    spostorder, tree = makepostorder(tree)
    epostorder = Order(n)
    firstsupernodelist = Vector{Int}(undef, m)
    lastsupernodelist = Vector{Int}(undef, m) 

    i₂ = 0

    for j in 1:m
        supernode = supernodelist[spostorder[j]]
        i₁ = i₂ + 1
        i₂ = i₂ + length(supernode)
        firstsupernodelist[j] = i₁
        lastsupernodelist[j] = i₂ 
        epostorder[i₁:i₂] .= supernode
    end 

    order = Order(order[epostorder], epostorder.index[order.index])
    firstancestorlist = [epostorder.index[firstancestorlist[spostorder[1:m - 1]]]; 0]
    subtreelist = spostorder.index[subtreelist[epostorder]]
    width = maximum(outdegreelist)

    EliminationTree(
        order,
        tree,
        firstancestorlist,
        firstsupernodelist,
        lastsupernodelist,
        subtreelist,
        width) 
end


# Construct a supernodal elimination tree, first computing an elimination order.
function EliminationTree(
    graph::AbstractGraph,
    algorithm::EliminationAlgorithm=DEFAULT_ELIMINATION_ALGORITHM,
    supernode::Supernode=DEFAULT_SUPERNODE)

    order = Order(graph, algorithm)
    EliminationTree(graph, order, supernode)
end


# Get the number of nodes in a supernodal elimination tree.
function Base.length(stree::EliminationTree)
    length(stree.tree)
end


# Get the width of a supernodal elimination tree.
function getwidth(stree::EliminationTree)
    stree.width
end


# Get the supernode at node i.
function getsupernode(stree::EliminationTree, i::Integer)
    i₁ = stree.firstsupernodelist[i]
    i₂ = stree.lastsupernodelist[i]
    stree.order[i₁:i₂]
end


# Get the highest node containing a vertex v.
function getsubtree(stree::EliminationTree, v::Integer)
    stree.subtreelist[stree.order.index[v]]
end


# Get the highest node containing vertices vs.
function getsubtree(stree::EliminationTree, vs::AbstractVector)
    stree.subtreelist[minimum(stree.order.index[vs])]
end


# Get the level of node i.
function getlevel(stree::EliminationTree, i::Integer)
    getlevel(stree.tree, i)
end


# Evaluate whether node i₁ is a descendant of node i₂.
function AbstractTrees.isdescendant(stree::EliminationTree, i₁::Integer, i₂::Integer)
    isdescendant(stree.tree, i₁, i₂)
end


# Draw a supernodal elimination tree.
function GraphPlot.gplot(
    stree::EliminationTree;
    nodelabel::Bool=false,
    nodefillc=RGB(0.1,0.1,0.1),
    rootfillc=RGB(0.882911,0.359638,0.360092),
    kwargs...)

    nl = nodelabel
    nc = parse(RGB{Float64}, nodefillc)
    rc = parse(RGB{Float64}, rootfillc)

    n = length(stree)
    graph = Graph(n)
    nodelabel = Vector{Vector{Int}}(undef, n)
    nodefillc = Vector{RGB{Float64}}(undef, n)

    nodelabel[n] = getsupernode(stree, n)
    nodefillc[n] = rc 

    for i in n - 1:-1:1
        add_edge!(graph, i, parentindex(stree, i))
        nodelabel[i] = getsupernode(stree, i)
        nodefillc[i] = nc
    end

    nodelabel = nl ? nodelabel : nothing
    nodelabelc = "white"

    gplot(graph; nodelabel, nodefillc, nodelabelc, kwargs...)
end


# Compute the supernodes, parent function, and first ancestor of a
# supernodal elimination tree.
#
# The complexity is
# 𝒪(n)
# where n = |V|.
#
# doi:10.1561/2400000006
# Algorithm 4.1: Maximal supernodes and supernodal elimination tree.
function makestree(etree::Tree, outdegrees::AbstractVector, supernode::Supernode)
    n = length(etree)
    sbt = Vector{Int}(undef, n)
    snd = Vector{Int}[]
    q = Int[]
    a = Int[]

    for v in 1:n
        w′ = findchild(etree, outdegrees, v, supernode)
        
        if isnothing(w′)
            i = length(snd) + 1
            sbt[v] = i
            push!(snd, [v])
            push!(q, 0)
            push!(a, 0)
        else
            i = sbt[w′]
            sbt[v] = i
            push!(snd[i], v)
        end

        for w in childindices(etree, v)
            if w !== w′
                j = sbt[w]
                q[j] = i
                a[j] = v
            end
        end
    end

    snd, sbt, q, a
end


# Find a child w of v such that
# v ∈ snd(w).
# If no such child exists, return nothing.
function findchild(etree::Tree, outdegrees::AbstractVector, v::Integer, ::Supernode) end


function findchild(etree::Tree, outdegrees::AbstractVector, v::Integer, ::MaximalSupernode)
    for w in childindices(etree, v)
        if outdegrees[w] == outdegrees[v] + 1
            return w
        end
    end
end


function findchild(etree::Tree, outdegrees::AbstractVector, v::Integer, ::FundamentalSupernode)
    ws = childindices(etree, v)

    if length(ws) == 1
        w = only(ws)

        if outdegrees[w] == outdegrees[v] + 1
            return w
        end
    end
end


# Compute the row and column counts of a graph's elimination graph.
#
# The complexity is
# 𝒪(m α(m, n))
# where m = |E|, n = |V|, and α is the inverse Ackermann function.
#
# doi:10.1137/S089547989223692
# Figure 3: Implementation of algorithm to compute row and column counts
function getdegrees(graph::AbstractGraph, order::Order, etree::Tree)
    n = nv(graph)
    forest = IntDisjointSets(n)   
    rvert = Vector{Int}(undef, n)
    index = Vector{Int}(undef, n)
    rvert .= index .= 1:n

    function FIND(p)
        index[find_root!(forest, p)]
    end

    function UNION(u, v)
        w = max(u, v)
        rvert[w] = root_union!(forest, rvert[u], rvert[v])
        index[rvert[w]] = w
    end
    
    postorder, etree = makepostorder(etree)
    graph = DiGraph(graph, compose(postorder, order))
    prev_p = Vector{Int}(undef, n)
    prev_nbr = Vector{Int}(undef, n)
    rc = Vector{Int}(undef, n)
    wt = Vector{Int}(undef, n)

    for u in 1:n
        prev_p[u] = 0
        prev_nbr[u] = 0
        rc[u] = 1
        wt[u] = isempty(childindices(etree, u))
    end

    for p in 1:n
        if p != n
            wt[parentindex(etree, p)] -= 1
        end

        for u in neighbors(graph, p)
            if getfirstdescendant(etree, p) > prev_nbr[u]
                wt[p] += 1
                p′ = prev_p[u]
                
                if p′ == 0
                    rc[u] += getlevel(etree, p) - getlevel(etree, u)
                else
                    q = FIND(p′)
                    rc[u] += getlevel(etree, p) - getlevel(etree, q)
                    wt[q] -= 1
                end
    
                prev_p[u] = p
            end

            prev_nbr[u] = p
        end

        if p != n
            UNION(p, parentindex(etree, p))
        end
    end

    cc = wt

    for v in 1:n - 1
        cc[parentindex(etree, v)] += cc[v]
    end

    indegrees = rc[postorder.index] .- 1
    outdegrees = cc[postorder.index] .- 1
    indegrees, outdegrees
end


##########################
# Indexed Tree Interface #
##########################


function AbstractTrees.rootindex(stree::EliminationTree)
    rootindex(stree.tree)
end


function AbstractTrees.parentindex(stree::EliminationTree, i::Integer)
    parentindex(stree.tree, i)
end


function AbstractTrees.childindices(stree::EliminationTree, i::Integer)
    childindices(stree.tree, i)
end


function AbstractTrees.NodeType(::Type{IndexNode{EliminationTree, Int}})
    HasNodeType()
end


function AbstractTrees.nodetype(::Type{IndexNode{EliminationTree, Int}})
    IndexNode{EliminationTree, Int}
end
