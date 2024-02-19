# A total ordering of the numbers {1, ..., n}.
struct Order <: AbstractVector{Int}
    order::Vector{Int}
    index::Vector{Int}
end


# Given a vector σ, construct the order ≺, where
#   σ(i₁) ≺ σ(i₂)
# if
#   i₁ < i₂.
function Order(order::AbstractVector)
    n = length(order)
    index = Vector{Int}(undef, n)

    for i in 1:n
        index[order[i]] = i
    end    

    Order(order, index)
end


# Construct an empty order of length n.
function Order(n::Integer)
    order = Vector{Int}(undef, n)
    index = Vector{Int}(undef, n)
    Order(order, index)
end


# Construct an elimination order using the reverse Cuthill-McKee algorithm. Uses
# CuthillMcKee.jl.
function Order(graph::AbstractGraph, ::CuthillMcKeeJL_RCM)
    order = CuthillMcKee.symrcm(adjacency_matrix(graph))
    Order(order)
end


# Construct an elimination order using the approximate minimum degree algorithm. Uses
# AMD.jl.
function Order(graph::AbstractGraph, ::AMDJL_AMD)
    order = AMD.symamd(adjacency_matrix(graph))
    Order(order)
end


# Construct an elimination order using the nested dissection heuristic. Uses Metis.jl.
function Order(graph::AbstractGraph, ::MetisJL_ND)
    order, index = Metis.permutation(graph)
    Order(order, index)
end


# Construct an elimination order using the maximum cardinality search algorithm.
function Order(graph::AbstractGraph, ::MCS)
    order, index = mcs(graph)
    Order(order, index)
end


function compose(order₁::Order, order₂::Order)
    order = order₂.order[order₁.order]
    index = order₁.index[order₂.index]
    Order(order, index)
end


# Evaluate whether
# n₁ < n₂
# in the given order.
function Base.isless(order::Order, n₁::Integer, n₂::Integer)
    order.index[n₁] < order.index[n₂]
end


# Compute a vertex elimination order using the maximum cardinality search algorithm.
#
# The complexity is
# 𝒪(m + n),
# where m = |E| and n = |V|. 
#
# https://doi.org/10.1137/0213035
# Maximum cardinality search
function mcs(graph::AbstractGraph)
    n = nv(graph)
    α = Vector{Int}(undef, n)
    α⁻¹ = Vector{Int}(undef, n)
    size = Vector{Int}(undef, n)
    set = Vector{Vector{Int}}(undef, n)

    set .= [[]]
    size .= 1
    append!(set[1], vertices(graph))

    i = n
    j = 1

    while i >= 1
        v = pop!(set[j])
        α[v] = i
        α⁻¹[i] = v
        size[v] = 0

        for w in neighbors(graph, v)
            if size[w] >= 1
                deletesorted!(set[size[w]], w)
                size[w] += 1
                insertsorted!(set[size[w]], w)
            end
        end

        i -= 1
        j += 1

        while j >= 1 && isempty(set[j])
            j -= 1
        end
    end

    α⁻¹, α
end


############################
# AbstractVector Interface #
############################


function Base.size(order::Order)
    (length(order.order),)
end


function Base.getindex(order::Order, i::Integer)
    order.order[i]
end


function Base.setindex!(order::Order, v::Integer, i::Integer)
    order.order[i] = v
    order.index[v] = i
    v
end


function Base.IndexStyle(::Type{Order})
    IndexLinear()
end
