###################################################################
# CategoricalTensorNetworks.jl                                    #
# https://github.com/AlgebraicJulia/CategoricalTensorNetworks.jl/ #
###################################################################


function contract_tensor_network(d::UndirectedWiringDiagram,
                                 tensors::AbstractVector{<:AbstractArray})
    @assert nboxes(d) == length(tensors)
    juncs = [junction(d, ports(d, b)) for b in boxes(d)]
    j_out = junction(d, ports(d, outer=true), outer=true)
    contract_tensor_network(tensors, juncs, j_out)
end


function contract_tensor_network(tensors::AbstractVector{<:AbstractArray{T}},
                                 juncs::AbstractVector, j_out) where T
    # Handle important binary case with specialized code.
    if length(tensors) == 2 && length(juncs) == 2
        return contract_tensor_network(Tuple(tensors), Tuple(juncs), j_out)
    end

    jsizes = Tuple(infer_junction_sizes(tensors, juncs, j_out))
    juncs, j_out = map(Tuple, juncs), Tuple(j_out)
    C = zeros(T, Tuple(jsizes[j] for j in j_out))
    for index in CartesianIndices(jsizes)
        x = one(T)
        for (A, junc) in zip(tensors, juncs)
            x *= A[(index[j] for j in junc)...]
        end
        C[(index[j] for j in j_out)...] += x
    end
    return C
end


function contract_tensor_network( # Binary case.
    (A, B)::Tuple{<:AbstractArray{T},<:AbstractArray{T}},
    (jA, jB), j_out) where T
    jsizes = Tuple(infer_junction_sizes((A, B), (jA, jB), j_out))
    jA, jB, j_out = Tuple(jA), Tuple(jB), Tuple(j_out)
    C = zeros(T, Tuple(jsizes[j] for j in j_out))
    for index in CartesianIndices(jsizes)
        C[(index[j] for j in j_out)...] +=
            A[(index[j] for j in jA)...] * B[(index[j] for j in jB)...]
    end
    return C
end


function infer_junction_sizes(tensors, juncs, j_out)
    @assert length(tensors) == length(juncs)
    njunc = maximum(Iterators.flatten((Iterators.flatten(juncs), j_out)))
    jsizes = fill(-1, njunc)
    for (A, junc) in zip(tensors, juncs)
        for (i, j) in enumerate(junc)
            if jsizes[j] == -1
                jsizes[j] = size(A, i)
            else
                @assert jsizes[j] == size(A, i)
            end
        end
    end
    @assert all(s >= 0 for s in jsizes)
    jsizes
end

