

function inner_blas2(x::Array{Float32,1}, y::Array{Float32,1}, dist::Array{Float32,1}, len::UInt64)
    @fastmath @inbounds for i = 1:len
        #println("left = ", (triang(i-1)+1), "\tright = ", triang(i))
        fragment!(view(dist, (triang(i-1)+1):triang(i)), view(x, 1:i), view(y, 1:i), x[i], y[i])
    end
end


function blas_stat_st2(indx, step::Float64, maxrng)
    F = Tuple.(indx)
    x, y= Float32.(first.(F)), Float32.(last.(F))
    dist = Vector{Float32}(undef,triang(length(x)))
    len::UInt64 = length(x)
    inner_blas2(x, y, dist, len)
    #print(length(fit(Histogram, dist, 0:step:(maxrng)).weights))
    return fit(Histogram, dist, 0:step:(maxrng)).weights
end

function fragment!(dist, arrx, arry, xi, yi)
    dist .= sqrt.((arrx .- xi) .^2 .+ (arry .- yi) .^2)
end

q = rand(Bool, 200,100)

v = findall(x->x==true, q)

function cart2arr(v)
    x, y = Float32.(first.(Tuple.(v))), Float32.(last.(Tuple.(v)))
end
function test(v)
    for i in 1:10
        x, y = Float32.(first.(Tuple.(v))), Float32.(last.(Tuple.(v)))
    end
end


#@btime ReconTPCF.blas_stat_st2(v, 1.0, 100)
x,y = cart2arr(v)

m = 3
n = 8
nx = 10000
ny = 8

a = rand(Float64, m, nx)
b = rand(Float64, m, ny)

#@btime fit(Histogram, vec(pairwise(Euclidean(), a, dims=2)))


using NearestNeighbors

using ReconTPCF

data = rand(1:1.0:100, 2, 10^4)
kdtree = KDTree(data; leafsize = 4)
NearestNeighbors.show(kdtree)

kdtree.hyper_rec
