
function fragment!(dist, arrx, arry, xi, yi)
    dist .= sqrt.((arrx .- xi) .^2 .+ (arry .- yi) .^2)
end

function inner_blas2(x::Array{Float32,1}, y::Array{Float32,1}, dist::Array{Float32,1}, len::UInt64)
    @fastmath @inbounds for i = 1:len
        #println("left = ", (triang(i-1)+1), "\tright = ", triang(i))
        fragment!(view(dist, (triang(i-1)+1):triang(i)), view(x, 1:i), view(y, 1:i), x[i], y[i])
    end
end

"""
    blas_stat_st2(indx, step, maxrng)

For a set of N points, compute the unique distance between all possible pairs of
points. Bins distance measurements with a histogram of bin width step and length maxrng
Single threaded implementation, greedy with memory. Good for smaller computations
"""
function blas_stat_st2(indx, step::Float64, maxrng)
    F = Tuple.(indx)
    x, y= Float32.(first.(F)), Float32.(last.(F))
    dist = Vector{Float32}(undef,triang(length(x)))
    len::UInt64 = length(x)
    inner_blas2(x, y, dist, len)
    #print(length(fit(Histogram, dist, 0:step:(maxrng)).weights))
    return fit(Histogram, dist, 0:step:(maxrng)).weights
end

using OnlineStats
using LinearAlgebra
using StatsBase

a = rand(N)
b = zeros(N, N)
b .= a .- a'
fit(Histogram, vec(LowerTriangular(b)), 0:0.1:1).weights

function compute_bbox(data::AbstractVector{V}) where {V <: AbstractVector}
    T = eltype(V)
    n_dim = length(V)
    maxes = Vector{T}(undef, n_dim)
    mins = Vector{T}(undef, n_dim)
    @inbounds for j in 1:length(V)
        dim_max = typemin(T)
        dim_min = typemax(T)
        for k in 1:length(data)
            dim_max = max(data[k][j], dim_max)
            dim_min = min(data[k][j], dim_min)
        end
        maxes[j] = dim_max
        mins[j] = dim_min
    end
    return (mins, maxes)
end


using StaticArrays, Tullio, LoopVectorization, BenchmarkTools, LinearAlgebra, StatsBase

#s = CartesianIndices(zeros(100,100))


s = findall(rand(Bool, 200, 100))
@btime ReconTPCF.blas_stat5(s, 100)
f = Tuple.(s)
q = reinterpret(SVector{2,Int}, s)
@btime compute_bbox(q)
@btime norm.(q)

const arr = rand(1:1.0:100, 2, 10000)
out = zeros(10^4, 10^4)
@btime @tullio $out[i,j] = sqrt(($arr[1,i] - $arr[1,j])^2 + ($arr[2,i] - $arr[2,j])^2) #This business is very fast
@tullio out[i,j] = sqrt((arr[1,i] - arr[1,j])^2 + (arr[2,i] - arr[2,j])^2)
w = copy(vec(LowerTriangular(out)))
@btime
function parhist(a, bins)
#let
    #a = Threads.@spawn fit(Histogram, vec(LowerTriangular(out))[], 0:1.0:100).weights
#end
    nt = Threads.nthreads()
    len = length(a)
    seg = div(len, nt)
    tasks = [Threads.@spawn fit(Histogram, view(a, ((i-1)*seg)+1:i*seg), bins).weights for i in 1:nt-1]
    out = [fetch(t) for t in tasks]
end
using CUDA
bins = 0:1.0:100
@benchmark parhist(w, bins)
@btime fit(Histogram, vec(LowerTriangular(out)), bins)
function kernel4(x)
    tid = threadIdx().x
    shared = @cuStaticSharedMem(Float32, 4)
    fill!(shared, 1f0)
    sync_threads()
    CUDA.atomic_add!(pointer(shared, tid), shared[tid + 2])
    sync_threads()
    return
end


using CUDAnative
using CuArrays

# CPU
function hist_cpu!(hist, δ, idx)
    Threads.@threads for j in 1:size(idx,2)
        @inbounds for i in 1:size(idx,1)
            hist[idx[i], j] += δ[i,j]
        end
    end
    return
end

# GPU
function hist_gpu!(h::CuMatrix{T}, x::CuArray{T}, id::CuArray{Int}; MAX_THREADS=256) where {T<:AbstractFloat}
    function kernel!(h::CuDeviceArray{T}, x::CuDeviceArray{T}, id)
        i = threadIdx().x + (blockIdx().x - 1) * blockDim().x
        j = threadIdx().y + (blockIdx().y - 1) * blockDim().y
        @inbounds if i <= size(id, 1) && j <= size(h, 2)
            k = Base._to_linear_index(h, id[i,j], j)
            CUDAnative.atomic_add!(pointer(h, k), x[i,j])
        end
        return
    end
    thread_i = min(MAX_THREADS, size(id, 1))
    thread_j = min(MAX_THREADS ÷ thread_i, size(h, 2))
    threads = (thread_i, thread_j)
    blocks = ceil.(Int, (size(id, 1), size(h, 2)) ./ threads)
    CuArrays.@cuda blocks=blocks threads=threads kernel!(h, x, id)
    return h
end

nbins = 20
ncol = 100
items = Int(1e6)
hist = zeros(Float32, nbins, ncol)
δ = rand(Float32, items, ncol)
idx = rand(1:nbins, items, ncol)

hist_gpu = CuArray(hist)
δ_gpu = CuArray(δ)
idx_gpu = CuArray(idx)

@time hist_cpu!(hist, δ, idx)
# 0.021154 seconds (64 allocations: 7.063 KiB)
@CuArrays.time hist_gpu!(hist_gpu, δ_gpu, idx_gpu, MAX_THREADS=1024)
# function get_bbox(indx)
#     F = Tuple.(indx)
#     x, y= Float32.(first.(F)), Float32.(last.(F))
#     minx = 0
#     miny = 0
#     maxx = 0
#     maxy = 0
#     for k in 1:length(x)
