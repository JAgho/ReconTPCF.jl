
using Distances
import Distances: Metric, result_type, eval_reduce, eval_end, eval_op, eval_start, evaluate, parameters
const MinkowskiMetric = Union{Euclidean,Chebyshev,Cityblock,Minkowski,WeightedEuclidean,WeightedCityblock,WeightedMinkowski}
struct KDNode{T}
    lo::T           # The low boundary for the hyper rect in this dimension
    hi::T           # The high boundary for the hyper rect in this dimension
    split_val::T    # The value the hyper rectangle was split at
    split_dim::Int  # The dimension the hyper rectangle was split at
end

struct KDTree{V <: AbstractVector,M <: MinkowskiMetric,T}
    data::Vector{V}
    indices::Vector{Int}
    metric::M
    nodes::Vector{KDNode{T}}
    reordered::Bool
end

a = KDTree()
