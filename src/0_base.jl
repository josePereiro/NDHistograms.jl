# TODO: integrate with https://juliastats.org/StatsBase.jl/stable/empirical/#StatsBase.NDHistogram
## ------------------------------------------------------------
# Work as dictionary, key are coordinates on the space, values are the counts. 
export NDHistogram
struct NDHistogram{elT}
    supports::Vector                # dimention support
    dim_names::Dict                 # dimention names <=> support index
    count_dict::Dict{elT, Int}      # space point => count
    function NDHistogram(ss::Pair...)
        isempty(ss) && error("Empty space definition.")
        supports = [last(si) for si in ss]
        dim_names = Dict()
        for (i, (name, _)) in enumerate(ss)
            dim_names[i] = name
            dim_names[name] = i
        end
        elT = Tuple{(eltype(last(si)) for si in ss)...}
        return new{elT}(supports, dim_names, Dict())
    end
end

# ------------------------------------------------------------
import Base.show
function Base.show(io::IO, h::NDHistogram{elT}) where elT
    println(io, "NDHistogram{", elT, "} with ", length(h), " non zero point(s)")
    for i in eachindex(h.supports)
        println(io, repr(h.dim_names[i]), " => ", repr(h.supports[i]))
    end
end

# ------------------------------------------------------------
import Base.haskey
Base.haskey(h::NDHistogram, v::Tuple) = haskey(h.count_dict, v)
import Base.getindex
Base.getindex(h::NDHistogram, v::Tuple) = getindex(h.count_dict, v)
import Base.get
Base.get(dflt::Function, h::NDHistogram, v::Tuple) = get(dflt, h.count_dict, v)
Base.get(h::NDHistogram, v::Tuple, dflt) = get(h.count_dict, v, dflt)
import Base.keys
Base.keys(h::NDHistogram) = keys(h.count_dict)
Base.keys(h::NDHistogram, dims) =
    (v[dimindex(h, dims)] for v in keys(h.count_dict))
import Base.values
Base.values(h::NDHistogram) = (w::Int for w in values(h.count_dict))

import Base.similar
function Base.similar(h0::NDHistogram, dims = Colon())
    dims = _col_dimindex(h0, dims)
    h1 = NDHistogram([
        dimname(h0, dim) => support(h0, dim)
        for dim in dims
    ]...)
    return h1
end

import Base.length
Base.length(h0::NDHistogram) = length(h0.count_dict)


## ------------------------------------------------------------
# name dimentions
export dimindex
dimindex(h::NDHistogram, idx::String) = h.dim_names[idx]
dimindex(h::NDHistogram, ::Colon) = eachindex(h.supports)
dimindex(::NDHistogram, idx) = idx
function _col_dimindex(h0::NDHistogram, dims)
    dims = dimindex(h0, dims)
    dims = dims isa Integer ? (dims:dims) : dims
    return dims
end

export dimname
dimname(h::NDHistogram) = [h.dim_names[i] for i in eachindex(h.supports)]
dimname(h::NDHistogram, dim) = dimname(h)[dimindex(h, dim)]

import Distributions.support
export support
Distributions.support(h::NDHistogram) = h.supports
Distributions.support(h::NDHistogram, dim) = support(h)[dimindex(h, dim)]

# for computing the discretization step
export delta_volume
delta_volume(r::AbstractRange) = step(r)
delta_volume(r) = mean(diff(r))
delta_volume(r::DataType) = one(r) # TODO: test this, think about it (it is the case where the support if all the space)

delta_volume(h::NDHistogram) = prod(delta_volume(s) for s in support(h))

## .-- .-. . .- .---. . ...- -- - --. ..- 
# assumes homogenous support
# TODO: implement non homeneous version
import Distributions.entropy
export entropy
function Distributions.entropy(h::NDHistogram)
    Σ_p_log_p = 0.0
    sum_w = sum(values(h))
    dx = delta_volume(h)
    for w in values(h)
        p = w / sum_w
        Σ_p_log_p += p * log(p / dx)
    end
    return -Σ_p_log_p 
end

# NOTE: assumes homogenoeus support
export pmf
function pmf(h::NDHistogram, v::Tuple; 
        sum_w = sum(values(h))
    )
    w = get(h, v, 0)
    return w / sum_w
end

# NOTE: assumes homogenoeus support
import Distributions.mean
function Distributions.mean(h::NDHistogram, dim)
    dim = dimindex(h, dim)
    sum_w = sum(values(h))
    Σvp = zeros(length(dim))
    for (v, w) in h.count_dict
        Σvp .+= v[dim] .* (w/sum_w)
    end
    return Tuple(Σvp)
end


# ------------------------------------------------------------
function _find_nearest(x::Real, x0::Real, dx::Real)
    iszero(dx) && return 1
    i, d = divrem(x - x0, dx)
    # @show i, d
    return d < (dx / 2) ? Int(i)+1 : Int(i)+2
end

function _find_nearest(x::Real, r::AbstractRange)
    x0, x1 = extrema(r)
    x0 > x && return firstindex(r)
    x1 < x && return lastindex(r)
    return _find_nearest(x, x0, step(r))
end

_sqerr(ri, x) = (x - ri)^2
_sqerr(x) = Base.Fix2(_sqerr, x)

_find_nearest(x::Real, r::AbstractVector) = last(findmin(_sqerr(x), r))

# ------------------------------------------------------------
# maps from Space -> bin
# eltype(Space) must return v's type
descretize(::Type{T}, v::T) where {T<:Number} = T(v)
descretize(S::AbstractRange, v) = getindex(S, _find_nearest(v, S))
descretize(S::AbstractVector, v) = getindex(S, _find_nearest(v, S))
descretize(F::Function, v) = F(v) # custom mapping
descretize(::T, n::T) where T = n # identity

# ------------------------------------------------------------
import Base.count!
function count!(h::NDHistogram, v::Tuple, w = 1)
    v = tuple([
        descretize(S, vi) for (S, vi) 
        in zip(h.supports, v)
    ]...)
    if haskey(h.count_dict, v)
        h.count_dict[v] += w
    else
        h.count_dict[v] = w
    end
    return h
end

const HISTOGRAMS_LK = ReentrantLock()
function lk_count!(h::NDHistogram, v::Tuple, w = 1)
    lock(HISTOGRAMS_LK) do
        count!(h, v, w)
    end
end

# ------------------------------------------------------------
# Merge histograms
import Base.merge!
function Base.merge!(h0::NDHistogram, h1::NDHistogram, hs::NDHistogram...)
    count_dict0 = h0.count_dict
    for (x, c) in h1.count_dict
        get!(count_dict0, x, 0)
        count_dict0[x] += c
    end
    for hi in hs
        for (x, c) in hi.count_dict
            get!(count_dict0, x, 0)
            count_dict0[x] += c
        end
    end
    return h0
end

## ------------------------------------------------------------
# TODO: find a better name
export rebin
function rebin(h0::NDHistogram, ss1::Pair...)
    ss1 = Dict(ss1...)
    for (n0, s0) in zip(dimname(h0), support(h0))
        haskey(ss1, n0) && continue
        ss1[n0] = s0
    end
    h1 = NDHistogram(pairs(ss1)...)
    # recount!
    for (v, w) in zip(keys(h0), values(h0))
        count!(h1, v, w)
    end
    return h1
end

## ------------------------------------------------------------
export marginal
function marginal(h0::NDHistogram, dims)
    dims = _col_dimindex(h0, dims)
    h1 = similar(h0, dims)
    for (v, w) in zip(keys(h0, dims), values(h0))
        count!(h1, v, w)
    end
    return h1
end

# ------------------------------------------------------------
import Base.filter
function Base.filter(f::Function, h0::NDHistogram)
    h1 = similar(h0)
    for (v, w) in h0.count_dict
        f(v) || continue
        count!(h1, v, w)
    end
    return h1
end

# ------------------------------------------------------------
# Give a sample vector with similar NDHistogram
# scale: scale back the sample vector size
# ex: [1,2] ~ [1,1,2,2] both vector has the same normilize histogram, 
# but rhs is half the side of lds
# It is great for plotting and reduce the number of points
# julia> samples = resample(h0, 2; scale)
# julia> lines!(ax, eachindex(samples) ./ scale, sort(samples))
export resample
function resample(h::NDHistogram, dims; scale = 1.0)
    _keys = keys(h, dims)
    _values = values(h)
    samples = Vector{typeof(first(_keys))}()
    for (x, w) in zip(_keys, _values)
        w = floor(Int, w * scale)
        push!(samples, fill(x, w)...)
    end
    return samples
end

## ------------------------------------------------------------
function _rand(rng::AbstractRNG, h::NDHistogram, Z; tries = 1e18)
    for _ in 1:tries # TODO: set a timeout
        v, w = rand(rng, h.count_dict)
        p = w/Z
        rand(rng) <= p && return v
    end
    error("Sampling failed!")
end

import Base.rand
function Base.rand(rng::AbstractRNG, h::NDHistogram; tries = 1e18)
    Z = sum(values(h))
    return _rand(rng, h, Z; tries)
end
Base.rand(h::NDHistogram; kwargs...) = rand(default_rng(), h; kwargs...)

function Base.rand(rng::AbstractRNG, h::NDHistogram{T}, n::Int; tries = 1e18) where T
    Z = sum(values(h))
    return [_rand(rng, h, Z; tries) for _ in 1:n]
end
Base.rand(h::NDHistogram, n::Int; kwargs...) = rand(default_rng(), h, n; kwargs...)

## ------------------------------------------------------------