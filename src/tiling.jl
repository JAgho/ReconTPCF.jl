


using BenchmarkTools
#using ReconTPCF
using StaticArrays

#floor(Int, 1/200*100)
#pwd(np.ndarray[np.double_t, ndim=2] r, np.ndarray[np.double_t, ndim=1] histo, double rmax, int nbins):
function dist1(r, histo, rmax, nbins)
#Calculate histogram of distances between particles belonging to the same structure (=set of coordinates).
    @inbounds @fastmath for i = 1:(size(r)[1]-1)
        for j = i+1:size(r)[1]
            d=0.
            @simd for k in 1:2
                tmp=r[i,k]-r[j,k]
                d+=tmp*tmp
            end
            d=sqrt(d)
            #d=sqrt(((r[i] - r[j])**2).sum())
            rIndex=ceil(Int64, d/rmax*nbins)
            if rIndex>=nbins
                #print("\n rIndex>=nbins\n")
                continue
            end
            histo[rIndex]+=1
        end
    end
    return histo
end

function dist2(r1, r2, histo, rmax::Int64, nbins::Int64)
#Calculate histogram of distances between particles belonging to different structures(=sets of coordinates).

    @inbounds @fastmath for i in 1:size(r1)[1]
        for j in 1:size(r2)[1]
            d::Float64 = hypot(r1[i,1]-r2[j,1], r1[i,2]-r2[j,2])

            #for k in 1:2
            #    tmp=r1[i,k]-r2[j,k]
            #    d+=tmp*tmp
            #end
            #d=sqrt(d)
            rIndex=ceil(Int64, d/(rmax*nbins))
            if rIndex>=nbins
                #print "\n rIndex>=nbins\n"
                continue
            end
            histo[rIndex]+=1
        end
    end
    return histo
end

function dist3(r1, r2, histo, rmax, nbins)
    @inbounds @fastmath for i in 1:size(r1)[1]
        for j in 1:size(r2)[1]
            d=0.

            tmp1=r1[i,1]-r2[j,1]
            tmp2=r1[i,2]-r2[j,2]
            d=sqrt(tmp1*tmp1 + tmp2*tmp2)
            rIndex=ceil(Int, d/rmax*nbins)
            if rIndex>=nbins
                #print "\n rIndex>=nbins\n"
                continue
            end
            histo[rIndex]+=1
        end
    end
    return histo
end

function dist4(r1, r2, histo, rmax, nbins)
    @inbounds @fastmath for i in 1:size(r1)[1]
        for j in 1:size(r2)[1]
            d=0.

            t::SVector{2,Int} = r1 - r2
            #tmp1=r1[i,1]-r2[j,1]
            #tmp2=r1[i,2]-r2[j,2]
            d=sqrt((t[1]*t[1]) + (t[2]*t[2]))
            rIndex=ceil(Int, d/rmax*nbins)
            if rIndex>=nbins
                #print "\n rIndex>=nbins\n"
                continue
            end
            histo[rIndex]+=1
        end
    end
    return histo
end

function dist5(r1, r2, histo, rmax, nbins, r1len::UInt, r2len::UInt)
    @inbounds @fastmath for i in 1:r1len
        for j in 1:r2len
            d=0.
            for k in 1:2
                tmp=r1[i,k]-r2[j,k]
                d+=tmp*tmp
            end
            d=sqrt(d)
            rIndex=ceil(Int, d/rmax*nbins)
            #if rIndex>=nbins
                #print "\n rIndex>=nbins\n"
            #    continue
            #end
            histo[rIndex]+=1
        end
    end
    return histo
end

function tile_dist(points::Array{Float64,2}, histo::Array{Int64,1}, rmax::Int64, nbins::Int64)
    bsize::Int64 = 100
    N = size(points)[1]
    histo .= 0
    if N <= bsize
        return dist1(points, histo, rmax, nbins)
    else
        histos = [zeros(Int64, length(histo)) for i in 1:Threads.nthreads()]
        #dump(histos[1])
        Threads.@threads for i in 0:(div(N, bsize)-1)
                            dist1(view(points, (i*bsize)+1:((i+1)*bsize), :),
                                 view(histos, Threads.threadid())[1],
                                 rmax,
                                 nbins)
            print("\nirange = ", (i*bsize)+1, ":", ((i+1)*bsize))
            #make thread with self_routine(myblocks)
        end
        @sync for il in 1:bsize:N+1
            for jl in il+bsize:bsize:N-bsize+1
                println("il: ", il, "\t\tir: ", il+bsize-1)
                println("jl: ", jl, "\t\tjr: ", jl+bsize-1)
                Threads.@spawn dist2(view(points, il:il+bsize-1, :),
                        view(points, jl:jl+bsize-1, :),
                        view(histos, Threads.threadid())[1],
                        rmax,
                        nbins)
            end
        end
        for i in 1:Threads.nthreads() histo .+= view(histos, i)[1] end
    end
    return histo
end

using FLoops
points = rand(1100, 2)
points2 = points.+1
hist = zeros(Int, 100)

tile_dist(points, hist, 5, 100)
@trace dist1(points, zeros(Int, 100), 100, 100)
@code_warntype dist2(points, points.+1, zeros(Int, 100), 100, 100)


sum(hist)
dist1(points, hist, 5, 100)
print(f)
view(points, 10:20, :)
points[10:20, :]

function staircase_looper(a, bsize)
    N = length(a)
    for (il, ir) in zip(1:bsize:N-1, bsize:bsize:N-1)
        for (jl, jr) in zip(ir+1:bsize:N-1, ir+bsize:bsize:N-1)
            println("il: ", il, "\t\tir: ", ir)
            println("jl: ", jl, "\t\tjr: ", jr)
        end
    end
end

staircase_looper(zeros(500), 100)




@btime dist3(points, points2, hist, 100, 100)
for i in 1:100
    dist1(points,  hist, 100, 100)
end
@btime dist5(points, points2, hist, 100, 100)





function solvequadratic(a, b, c)
    d = sqrt(b^2 - 4a*c)
    (-b - d) / 2a, (-b + d) / 2a
end

solvequadratic(-0.5, 999.5, -90000.0)






function bench_dist(N)
    points = rand(1:1.0:100, N, 2)
    points2 = points.+1
    hist = zeros(Int, 100)
    return @benchmark dist2($points, $points2, $hist, 100, 100)
end
a = []
range = 1000:1000:10000
for i in range
    #print("for i ops = ", i*i/2 , " ")
    push!(a, bench_dist(i))
end
for (bench, i) in zip(a, range)
    println("for ", abs2(i)/2, " operations, ", minimum(bench).time/(abs2(i)/2))
end

#q = findall(ones(Bool, 10,10))
#s1 = reinterpret(SVector{2,Int}, q)
#t = zeros(Bool, 20,10)
#t[11:end,:] .= true
#z = reinterpret(SVector{2,Int}, findall(t))
#s2 = reinterpret(SVector{2,Int}, q)
#dist4(s1, s2, zeros(Int, 100), 100, 100)
#@inline function euclid(a, b) (a.-b)


@btime ReconTPCF.blas_stat5(q, 100)
