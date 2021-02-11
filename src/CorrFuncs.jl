using Images, FileIO
using StatsBase
using AverageShiftedHistograms
using DSP
using FFTW

#using CSV
#using DataFrames

#Pre processes an image
#takes a 2D array
#returns a 2D array


"""
    pre_proc(im)

convolve image with a disk-shaped structure element
"""
function pre_proc(im)
    img = Gray.(im) .> 0.05
    bw = Float64.(img)
    strel = [0 1 1 1 0;
             1 1 1 1 1;
             1 1 1 1 1;
             1 1 1 1 1;
             0 1 1 1 0]
    R = findall(x->x==true, strel)
    LocalFilters::erode!(bw, R)
    LocalFilters::dilate!(bw, R)
    #localfilter!(bw, :, min, 3)
    #localfilter!(bw, :, max, 3)
    im2 = nonzero(bw)
    return im2
end


"""
    nonzero(image)

evaluate where an image is
nonzero
"""
function nonzero(image)
    dims = size(image)
     im = Array{Bool}(undef, dims[1], dims[2])
    for i = 1:dims[1]
        for j = 1:dims[2]
            im[i,j] = !(image[i,j]>0)
        end
    end
    return im
end

#Sloppy euclidian distance computation
#takes a 1D array of CartesianIndex tuples
#returns a 1D array of Float32
function pdist(indices)

    dims = size(indices)[1]

    lens = Array{Float32}(undef, Nlen(size(indices)[1]+1))
    idx = 1

    for i = 1:dims[1]
            for j = i:dims[1]
                lens[idx] = cart(indices[i], indices[j])
                idx += 1
            end
    end

    return lens
end


#Computes the euclidian distance between a pair of CartesianIndex tuples
#takes a tuple of two CartesianIndex tuples
#returns a float
function cart(a, b)
    return sqrt((a[1]-b[1])^2 + (a[2]-b[2])^2)
end

#Computes distances between a list of indices, and histograms these
#values across a 1:100 interval with bins of width 1

#takes a list of CartesianIndices
#returns an AverageShiftedHistograms ash object (see methods on github)
#TODO make the range adjustable, add optimisation to shortcut for cluster functions
#(computing this many histograms is wasteful vs one big one)
function im_stat(indx)
    dims = size(indx)[1]
    lens = zeros(Float32, dims)
    primer = Array{Float32}([1])
    s1 = ash(primer, rng=0:1:100.)
    idx = 1
    @inbounds @fastmath  Threads.@threads for i = 1:dims[1]
        #a = indx[i]
            @inbounds @fastmath for j = i:((dims)[1])
                lens[j] = cart(indx[i], indx[j])
               #print("\n i = ", i, " j = ", j)
            end
        #idx += 1

        merge!(s1, ash(lens[idx:end], rng=0:1:100.))
    end

    return s1
end

function im_stat1(indx)
    dims = size(indx)[1]
    lens = zeros(Float32, dims)
    primer = Array{Float32}([1])
    s1 = ash(primer, rng=1:1:100.)
    idx = 1
    @inbounds @fastmath for i = 1:dims[1]
        #a = indx[i]
            @simd for j = idx:((dims)[1])
                lens[j] .= cart(indx[i], indx[j])
               #print("\n i = ", i, " j = ", j)
            end
        idx += 1

        merge!(s1, ash(lens[idx:end], rng=1:1:100.))
    end

    return s1
end

#Computes distances between a list of indices, and histograms these
#values across a 1:100 interval with bins of width 1
#does this after computing the 8-connectivity set of clusters
#does this for all clusters. Also computes the SN of the image if not provided

#takes an image
#returns an AverageShiftedHistograms ash object (see methods on github)
#TODO tuple this with SN as output
function tpcof(im, SN=0...)

    if SN == 0
        full_arr = trues((1:size(im)[1], 1:size(im)[2]))

        ind = findall(x->x==true, full_arr) #CartesianIndices((1:size(im)[1], 1:size(im)[2]))
        #ind = CartesianIndices((size(im)[1], size(im)[2]))
        Nind = length(ind)
        SN = (2 .* im_stat(ind).counts) ./ length(im)
        print("\n")

    end

    indices = findall(x->x==true, im)
    BN = (2 .* im_stat(indices).counts) ./ length(indices)
    S2 = (length(indices)/length(im)) .* (BN ./ SN)
    #print(length(indices), "    ", length(im) , "\n")

    #S2 =  (length(indices)/length(im)).*((BN.counts)./ (2 .* SN.counts ./))
    S2[1] = length(indices)/length(im)

    return S2, length(indices), SN, 0:1:100.
end


#Computes distances between a list of indices, and histograms these
#values across a 1:100 interval with bins of width 1
#Also computes the SN of the image if not provided

#takes an image
#returns an AverageShiftedHistograms ash object (see methods on github)
#TODO tuple this with SN as output
function tpclf(im, SN=0...)

    label = label_components(im)
    comps = component_subscripts(label)[2:end]
    bbox = component_boxes(label)
    primer = Array{Float32}([1])
    s1 = ash(primer, rng=1:1:100.)

    if SN == 0
        ind = CartesianIndices((1:size(im)[1], 1:size(im)[2]))
        Nind = length(ind)
        SN = im_stat(ind)
    end

    for cluster in comps
        merge!(s1, im_stat(cluster))
    end
    #print(size(xy(s1))
    C2 = (Nind./length(im).*(xy(s1)[2]./ xy(SN)[2]))

    return C2
end


"""
Compute a binary disk-shaped structure element
given radius r on the basis of (areal) interpolation

"""
function disk_strel(r)
    dim = (2*r)-1
    arr = zeros(Bool, dim, dim)
    coords = CartesianIndices(arr)
    #dump(coords)
    middle = (length(coords) ÷ 2) +1
    mid = coords[middle]
    for idx in coords
        if cart(idx, mid) < 5
            arr[idx] = true
        end

    end
    return arr
    #centre = CartesianIndices(arr[r,r])


end

function polydisp_iso_circs(dims, r, N)
    x = dims[1]
    y = dims[2]
    off = (r-1)
    im = zeros(Bool, x+(off*2), y+(off*2))
    strel = disk_strel(r)
    #carts = findall(x->x==true, view(im, off:x, off:y))
    for i= 1:N
        x1 = rand(off+1:x)
        y1 = rand(off+1:y)
        im[(x1-off):(x1+off), (y1-off):(y1+off)] .= im[(x1-off):(x1+off), (y1-off):(y1+off)] .| strel
        #print(x1,"\t", y1, "\n")
    end
    #dump(carts)
    #strel = disk_strel(r)
    #centres = sample(carts, 80)
    #for coord in centres
        #x = coord[1]
        #y = coord[2]
        #print(coord)
        #im[(x-off):(x+off), (y-off):(y+off)] = strel


return im
end


"""
Generate isotropically distributed circles that
may overlap. Keep the implied phase fraction
ϕ below 0.3 if possible

"""
function monodisp_circ(dims, R, N)
    x = dims[1]
    y = dims[2]
    im = zeros(Bool, dims)
    xrand = rand(1:x, N)
    yrand = rand(1:y, N)
    coords = hcat(xrand, yrand)
    for coord in zip(xrand, yrand)
        #print(coord, "\n")
        for i in 1:x, j in 1:y

            if cart((i,j), coord) < R
                im[i,j] = true
            end
        end
    end
    return im

end

function S2_analyt(r, R, phi)
    #v2dv1=2*(r>2*R)+2/pi*(pi+r/2/R.*sqrt(1-r.^2/4/R^2)-acos(r/2/R)).*(2*R>=r);
    squarebit = real(sqrt.(Complex.((1 .- r .^2 ./ 4 ./ (R^2)))))
    f = (r ./ 2 ./ R)
    cosbit = acos.(f .* (f .<= 1 )) .* (f .<= 1 )
    leftbit = 2 .* (r .> (2*R))
    rightbit = 2 .* R .>=r
    midbit = r ./ 2 ./ R
    total = leftbit .+ 2 ./ pi .* (pi .+ (midbit .* squarebit) .-  cosbit) .* rightbit

    eta=-log(phi);

    S2=exp.(-total .* eta);

    return S2
end


function S2_fft(I)
    n,m = size(I)
    I_fft = fft(I)
    abs.(fftshift(ifft(I_fft .* conj(I_fft)))) ./ (n*m)
end

function S2_fft_test1(I,P)
    n,m = size(I)
    I_fft = fft(I, P)
    abs.(fftshift(ifft(I_fft .* conj(I_fft)))) ./ (n*m)
end

function S2_fft_test(I)

    I_fft = rfft(I)

    n,m = size(I_fft)
    abs.(fftshift(ifft(I_fft .* conj(I_fft)))) ./ (n*m)
end

function C2_pdist(im, SN=0...)

    label = label_components(im)
    comps = component_subscripts(label)[2:end]
    bbox = component_boxes(label)
    primer = Array{Float32}([1])
    s1 = ash(primer, rng=0:1:100.)

    if SN == 0
        ind = CartesianIndices((size(im)))
        Nind = length(ind)
        SN = im_stat(ind)
    end

    for cluster in comps
        merge!(s1, blas_stat(cluster))
    end
    #print(size(xy(s1))
    C2 = (Nind./length(im).*(xy(s1)[2]./ xy(SN)[2]))
    for i= 1:length(C2)
    if isnan(C2[i])
        C2[i] = 0
    end
end

    return C2
end

function C2_fft_MT(I)

    l = label_components(I)
    bbox = component_boxes(l)[2:end]
    C2= zeros(size(I)[1], size(I)[2], Threads.nthreads())
    #A = zeros(Threads.nthreads())
    alloc = falses(size(I)[1], size(I)[2], Threads.nthreads())
    Threads.@threads for cluster in bbox
        #print(S2_fft(cluster))
        xbound = cluster[1][1]:cluster[2][1]
        ybound = cluster[1][2]:cluster[2][2]
        alloc[xbound, ybound, Threads.threadid()] = I[xbound, ybound]
        #print(size(alloc[:,:,Threads.threadid()]))
        #print("\n")
        #print(size(C2[:,:, Threads.threadid()]))
        C2[:,:, Threads.threadid()] = C2[:,:, Threads.threadid()] .+ S2_fft(alloc[:,:,Threads.threadid()])
        #if mod(i , 100) == 0
        #    imshow(C2)
        #end
        alloc[xbound, ybound, Threads.threadid()] .= false
        #C2 = C2 .+ S2_fft(cluster)
    end


    return C2
end

function C2_fft_MT_test(I)

    P = plan_fft(I)
    l = label_components(I)
    bbox = component_boxes(l)[2:end]
    C2= zeros(size(I)[1], size(I)[2], Threads.nthreads())
    alloc = falses(size(I)[1], size(I)[2], Threads.nthreads())

    Threads.@threads for cluster in bbox
        xbound = cluster[1][1]:cluster[2][1] #find bounds
        ybound = cluster[1][2]:cluster[2][2]
        alloc[xbound, ybound, Threads.threadid()] = I[xbound, ybound] #bind values to preallocated array
        C2[:,:, Threads.threadid()] = C2[:,:, Threads.threadid()] .+ S2_fft(alloc[:,:,Threads.threadid()])
        alloc[xbound, ybound, Threads.threadid()] .= false

    end


    return C2
end

function blas_stat(indx)
    F = Tuple.(indx)
    x, y= first.(F), last.(F)
    x1, y1 = copy(x), copy(y)
    len = size(x)[1]
    #print(len)
    xout = zeros(len)
    yout = zeros(len)
    dist = zeros(len)
    primer = Array{Float32}([1])
    s1 = ash(primer, rng=0:1:100.)
    @inbounds @fastmath Threads.@threads for i = 1:len
        dist .= sqrt.((x .- x1[1]) .^2 .+ (y .- y1[1]) .^2)
        ash!(s1, view(dist, i:len))
    end
    return s1
end

function C2_fft(I)
    l = label_components(I)
    bbox = component_boxes(l)[2:end]
    C2= zeros(size(I)[1], size(I)[2])
    #A = zeros(Threads.nthreads())
    alloc = falses(size(I)[1], size(I)[2])
    for cluster in bbox
        #print(S2_fft(cluster))
        xbound = cluster[1][1]:cluster[2][1]
        ybound = cluster[1][2]:cluster[2][2]
        alloc[xbound, ybound] = I[xbound, ybound]
        #print(size(alloc[:,:,Threads.threadid()]))
        #print("\n")
        #print(size(C2[:,:, Threads.threadid()]))
        C2 = C2 .+ S2_fft(alloc[:,:])
        #if mod(i , 100) == 0
        #    imshow(C2)
        #end
        alloc[xbound, ybound] .= false
        #C2 = C2 .+ S2_fft(cluster)
    end
    return C2
end


"""
    make_rand_im(philen, dims)

create a random binary array of size dims with
philen white pixels set within it. Also return
a coordinate list for white and black pixels
"""
function make_rand_im(philen, dims)
    idx = CartesianIndices((dims))
    #println(length(idx))
    wpix = sample(idx, philen, replace=false)
    #println(length(wpix))
    #println(wpix)
    bpix = setdiff(idx, wpix)
    im = zeros(Bool, (dims))

    for pix in wpix
        im[pix] = true
    end
    return im, (wpix, bpix)
end


"""
    get_rand_pix((pix))::Tuple{Array{CartesianIndex{2},1},Array{CartesianIndex{2},1}}

Select a random white and black pixel from a tupled
pair of lists of cartesian coordinates
"""
function get_rand_pix((pix))::Tuple{Array{CartesianIndex{2},1},Array{CartesianIndex{2},1}}
    wpick = [sample(pix[1])]
    bpick = [sample(pix[2])]
    #print(wpick, "\t", bpick, "\n")

    return wpick, bpick
end



function ediff_S2_fft(S2_guess,  S2)
    return sum((S2 .- S2_guess) .^2)
end

function ediff_C2(C2_guess, C2)
    return sum((((C2 .- C2_guess) .* 2) .^2))
end

function ediff_S2(S2_guess, S2)
    return sum((((S2 .- S2_guess) .* 2) .^2))
end

function ediff_test(C2_guess, C2, S2_guess,  S2)
    return (ediff_C2(C2_guess, C2), ediff_S2_fft(S2_guess, S2))
end

"""
    cluster_stat(clusters, step::Float64, maxrng)

Fast compution of S2 count for an array of arrays of cartesian indices.
Multithreaded in the sense that
"""
function cluster_stat(clusters, step::Float64, maxrng)
    s1_cache_threads = [zeros(Int64,maxrng) for i in 1:Threads.nthreads()];
    @inbounds @fastmath Threads.@threads for indx in clusters
        s1_cache_threads[Threads.threadid()] .+= blas_stat_st2(indx, step, maxrng)
    end
    f = zeros(Int64, length(s1_cache_threads[1]))
    for i in 1:Threads.nthreads(); f = f .+ s1_cache_threads[i]; end
    return f
end

function findlens(list, lower, upper)::Tuple{Int64,Int64}
    mymap = map( x->(lower<=length(x)<=upper), list)
    left = findfirst(isequal(1), mymap)
    right = findlast(isequal(1), mymap)
    #println("left = ", left, "\tright = ", right)
    if left === nothing
        return (0,0)
    else
    return (left,right)
    end
end

"""
    C2_pdist4(im, maxrng::Int64, SN::Array{Float64,1}=Array{Float64}(undef, 0))

Fast total computation of the C2 count for an array of arrays of cartesian indices
This version splits the problem by the type of cluster, selecting an appropriate
algorithm to suit the scale of each component.
"""
function C2_pdist4(im, maxrng::Int64, SN::Array{Float64,1}=Array{Float64}(undef, 0))

    comps = component_subscripts(label_components(im, trues(3,3)))[2:end]
    sort!(comps, by = x -> (length(x)))
    #println(length.(comps))
    singles = findlens(comps, 1, 1)
    #dump(singles)
    smalls = findlens(comps, 2, 20)
    bigs = findlens(comps, 21, length(im))
    #println("singles = ", singles, " smalls = ", smalls, " bigs = ", bigs)
    #bbox = component_boxes(label)
    #primer = Array{Float32}([])
    #maxrng = 100#min(100, size(im)[1], size(im)[2])
    BN = fit(Histogram, Array{Float32}([]), 0:1.:maxrng, closed=:left).weights

    if isempty(SN)#SN == 0.0
        SN = SN_comp(size(im),  maxrng)
        #ind = vec(CartesianIndices((size(im))))

        #SN = 2 .* blas_stat4(ind, 1.0, maxrng) ./ length(im)
    end
    #println("there are ", length(comps), " clusters")
    #println("SN = ", SN)
    #plot(SN)
    if singles !==(0,0) #Accumulate BN for all single-pixel clusters (dist=0 implicitly)
        BN[1] += 2*singles[2]
        #println("trig'ed\t", length(singles))
    end

    if smalls !==(0,0)
        BN .= (BN .+ (2 .* cluster_stat((comps[smalls[1]:smalls[2]]), 1.0, maxrng)))
    end

    if bigs!==(0,0)
        for cluster in (@view comps[bigs[1]:bigs[2]])
            BN .= (BN .+ (2 .* blas_stat4(cluster, 1.0, maxrng)))
            #break
        end
    end
    BN =  (BN ./ sum(im))
    C2 = ((sum(im)/length(im)) .* (BN ./ SN))
    for i= 1:length(C2)
        if isnan(C2[i])
            C2[i] = 0
        end
    end
    return C2, BN
end

"""
    SN_comp(dims::Tuple{Int64,Int64},  maxrng::Int64)

Wrapper for ```blas_stat5``` acting as a function barrier
"""
function SN_comp(dims::Tuple{Int64,Int64},  maxrng::Int64)::Array{Int64,1}
    2 .* blas_stat5(vec(CartesianIndices(dims)), maxrng)
end


"""
    inner_blas2(x, y, dist, len)

compute the lower triangular self-interaction matrix of a list of N items.
gives a flat vector of F32, of length (N(N-1)/2). Interaction here is computing
the L2 norm, but could be otherwise.
"""
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

function blas_stat_st2(x::Array{Float64,1}, y::Array{Float64,1}, step::Float64, maxrng)
    F = Tuple.(indx)
    x, y= Float32.(first.(F)), Float32.(last.(F))
    dist = Vector{Float32}(undef,triang(length(x)))
    len::UInt64 = length(x)
    inner_blas2(x, y, dist, len)
    #print(length(fit(Histogram, dist, 0:step:(maxrng)).weights))
    return fit(Histogram, dist, 0:step:(maxrng)).weights
end

function triang(n::UInt64)::UInt64
    return (n*(n+1))÷2
end

function triang(n::Int64)::Int64
    return (n*(n+1))÷2
end


"""
Compute the L2 norm between points (xi, yi) and the arrays (arrx, arry) and write
to preallocated array dist
"""
function fragment!(dist, arrx, arry, xi, yi)
    dist .= sqrt.((arrx .- xi) .^2 .+ (arry .- yi) .^2)
end


"""
    blas_stat4(indx, step::Float64, maxrng)::Array{Int64,1}

Computes a very fast multithreaded probabilistic S2 count. Employs kernel density
estimation in place of a true histogram. Gives hilarious results when used
for reconstruction.
"""
function blas_stat4(indx
    , step::Float64, maxrng)::Array{Int64,1}
    F = Tuple.(indx)
    x, y= Float32.(first.(F)), Float32.(last.(F))
    len = size(x)[1]
    primer = Array{Float32}([0])
    _s1_cache_threads = [ash(primer, rng=0:step:(maxrng-1), m =1) for i in 1:Threads.nthreads()];
    _dist_cache_threads = [Vector{Float32}(undef,len) for i in 1:Threads.nthreads()]

    @inbounds @fastmath Threads.@threads for i = 1:len
        fragment!(view(_dist_cache_threads, Threads.threadid())[1], x, y, x[i], y[i] )
        println(view(_dist_cache_threads[Threads.threadid()], i:len))
        ash!(view(_s1_cache_threads, Threads.threadid())[1], view(_dist_cache_threads[Threads.threadid()], i:len))
    end
    f = zeros(Int64, length(_s1_cache_threads[1].counts))
    for i in 1:Threads.nthreads(); f .= f .+ _s1_cache_threads[i].counts; end
    f[1] = f[1] - Threads.nthreads()
    return f
end


"""
Compute the S2 count (multithreaded, O(N) memory)

For a list of N CartesianIndex tuples, compute and bin the L2 norm between all pixel
pairs in the list.
"""
function blas_stat5(indx, maxrng)
    F = Tuple.(indx)
    x, y= Float32.(first.(F)), Float32.(last.(F))
    len = size(x)[1]
    primer = Array{Float32}([])
    cnt_thrd = [zeros(Int64,maxrng) for i in 1:Threads.nthreads()];
    dist_thrd = [zeros(Float32, len) for i in 1:Threads.nthreads()]

    @inbounds @fastmath Threads.@threads for i = 1:len
        fragment!(view(dist_thrd, Threads.threadid())[1], x, y, x[i], y[i] )
        view(cnt_thrd, Threads.threadid())[1] .+=  fit(Histogram, view(dist_thrd[Threads.threadid()], i:len), 0:1:maxrng).weights
    end
    f = zeros(Int64, length(cnt_thrd[1]))
    for i in 1:Threads.nthreads(); f = f .+ cnt_thrd[i]; end
    #f[1] = f[1] - Threads.nthreads()
    return f
end

"""
Get list of all surface pixels from a 2d binary image using convolution to discriminate
"""
function surf_opt2(im::Array{Bool,2})::Tuple{Array{CartesianIndex{2},1},Array{CartesianIndex{2},1}}
    edge_kern::Array{Int64} = [1 1 1; 1 10 1; 1 1 1]
    arrim::Array{Int64,2} = conv(im, edge_kern)
    return (findall(in(10:17), arrim), findall(in(1:8), arrim))
end

"""
Get list of all surface pixels from a 2d binary image using convolution to discriminate
Sanitised version to prevent the function returning nothing
"""
function surf_opt(im::Array{Bool,2})::Tuple{Array{CartesianIndex{2},1},Array{CartesianIndex{2},1}}
    edge_kern::Array{Int64} = [1 1 1; 1 10 1; 1 1 1]
    arrim::Array{Int64,2} = conv(im, edge_kern)[2:end-1, 2:end-1]
    #display(heatmap(arrim, color=:grays, aspect_ratio=1))
    #imshow(im)
    #imshow(arrim)
    bb::Array{CartesianIndex{2},1} = findall(in(1:8), arrim)
    wb::Array{CartesianIndex{2},1} = findall(in(10:17), arrim)

    if !isnothing(bb) & ! isnothing(wb)
        return (wb, bb)
    else
        return  (CartesianIndices(collect(1:1)), CartesianIndices(collect(1:1)))
    end
end

"""
Update tuple of 2 cartesian indices by swapping elements from each list at linear
index wpick (for 1) and bpick (for 0)
"""
function swap_pix(guess, pix, wpick, bpick)
    #println("swap_pix : \t\twpix length = ", length(pix[1])," bpix length = ", length(pix[2]) )
    guess[pix[1][wpick]] = false ##THIS IS THE DANGER ZONE
    guess[pix[2][bpick]] = true ########THIS IS NOT SYMMETRIC

    pix[1][wpick], pix[2][bpick] = pix[2][bpick], pix[1][wpick]

    return (pix)
end


"""
    find_equivalent(shortlist, pix, pick_s)

Computes equivalent position of unique entry in a shortlist to an entry in a master
list. Uses a much more sensible tupled pick - should use that everywhere
"""
function find_equivalent(shortlist, pix, pick_s)
    wcart = shortlist[1][pick_s[1]]
    bcart = shortlist[2][pick_s[2]]
    #println("find_equivalent : \twpick = ", pick_s[1], "\tbpick = ", pick_s[2])
    #println("find_equivalent : \twpix len = ", length(pix[1]), "\tbpix len = ", length(pix[2]))
    wrand = findfirst(x -> x == wcart, pix[1])
    brand = findfirst(x -> x == bcart, pix[2])
    #println("find_equivalent : \twrand = ", wrand, "\tbrand = ", brand)
    return (wrand, brand)
end

"""
Find the intersection between two tupled lists of cartesian indices
We use this to create a segment of a list that maintains ordering
"""
function surf_rand(surf_res, pix)
    return (intersect(surf_res[1], pix[1]), (intersect(surf_res[2], pix[2])))
end

"""
Selects random elements in a pair of lists
"""
function pick_rand_coord(list::Tuple{Array{CartesianIndex{2},1},Array{CartesianIndex{2},1}})::Tuple{Int64, Int64}
    wpick = rand(1:length(list[1]))
    #println("pick_rand_coord : \twlen = ", length(list[1]), "\tblen = ", length(list[2]))
    bpick = rand(1:length(list[2]))
    #println("pick_rand_coord : \twpick = ", wpick, "\tbpick = ", bpick)
    return wpick, bpick
end

"""
Fetch the region around a pixel selection; give a reduced window if near an edge
"""
function fetch_locale(im, pix, kern, wpick, bpick)
    re, be = size(im)[1], size(im)[2]
    wp = pix[1][wpick]
    bp = pix[2][bpick]

    #println(kern)
    #println("wp = ",wp)
    kern = check_edge(wp, re, be)
    #println(kern)
    #println(wp .+ kern)
    a = view(im, wp .+ kern)
    kern = check_edge(bp, re, be)
    #println("bp = ", bp)
    #println(kern)
    #println(bp .+ kern)
    b = view(im, bp .+ kern)
    return (a, b)
end


"""
Prevent fetch_locale from colliding with image edges by providing a reduced
window if an edge is collided with
"""
function check_edge(p, re, be)
    kern = CartesianIndices((-1:1, -1:1))
    l1 = -1
    r1 = 1
    l2 = -1
    r2 = 1
    if p[1] ==1
        l1 = 0
        #kern = kern .= CartesianIndices((0:1, -1:1))
    end
    if p[1] ==re
        r1 = 0
       # kern = CartesianIndices((-1:0, -1:1))
    end
    if p[2] ==1
        l2 = 0
        #kern = CartesianIndices((-1:1, 0:1))
    end
    if p[2] ==be
        r2 = 0
        #kern = CartesianIndices((-1:1, -1:0))
    end
    kern = CartesianIndices((l1:r1, l2:r2))


    return kern
end

"""
Find all unique values present in a small window into the region adjacent to a
pixel selection. Used to interrogate a cluster labelled image
"""
function get_region_names(regions)
    w1 = unique(regions[1][1])
    w2 = unique(regions[2][1])
    b1 = unique(regions[1][2])
    b2 = unique(regions[2][2])
    return ((w1, b1),(w2, b2))
end

"""
Deduct the contribution of a given set of pixels from a previously calculated
S2 count
"""
function subtract_cluster(hist, cluster)
    hist .= hist .- blas_stat_st2(cluster, 1.0, 100)
end

"""
Examine which clusters are present in multiple windows simultaneously to determine
which clusters are due to change from a pixel swap
"""
function find_unique_regions(names)
    olds = filter(x->x!=0,union(names[1][1], names[1][2]))
    news = filter(x->x!=0,union(names[2][1], names[2][2]))

    return (olds, news)
end

"""
Compute the C2 count in an image. Serves as a starting point for a reconstruction
"""
function C2_initialise(im, maxrng::Int64)
    comps = component_subscripts(label_components(im, trues(3,3)))[2:end]
    BN = zeros(Int64, maxrng)
    #for comp in comps
    #println(length(comp))
    #end
    BN .= (BN .+ 2 .* (cluster_stat((comps), 1.0, maxrng)))
    return BN
end

"""
Sloppy computation to find the contribution to C2 count from the old state and the
new, and hence the difference between them.
"""
function compute_change!(subtract, add, clusts_old, clusts_new, maxrng)
    for clust in clusts_old; subtract .= subtract .+ blas_stat5(clust, maxrng);  end
    for clust in clusts_new; add      .= add .+ blas_stat5(clust, maxrng); end
    return  add, subtract
end


"""
Eliminates NaN and Inf in an array
"""
function naninf(arr)
    for i=1:length(arr)
        if isnan(arr[i]) | isinf(arr[i])
            arr[i] = 0
        end
    end
    return arr
end


"""
Computes the C2 count change between a pair of images with one pixel swapped between them
"""
function update_C2_BN(BN, pix, original::Array{Int64, 2}, modified::Array{Int64, 2}, wpick, bpick, maxrng)
    kern = CartesianIndices((-1:1, -1:1))
    names = get_region_names((fetch_locale(original, pix, kern, wpick, bpick), fetch_locale(modified, pix, kern, wpick, bpick)))
    old, new = find_unique_regions(names)
    #println("old unique regions = ", old)
    #println("new unique regions = ", new)
    subtract = zeros(Int64, length(BN))
    add = zeros(Int64, length(BN))
    #println("BN is of length ",  length(BN))
    #clusts_old = component_subscripts(original) #we can replace this with findall for a particular index or indices of interest
    #clusts_new = component_subscripts(modified)
    clusts_old, clusts_new = get_clusters(old, new, original, modified)
    #println(size(clusts_new), "\t", size(clusts_old) ,"\t", clusts_new, "\t",clusts_old)
    add, subtract = compute_change!(subtract, add,  clusts_old, clusts_new, maxrng)

    return  (add .- subtract) .* 2
end

"""
Fuses clusters to reduce the scale of the problem
"""
function get_clusters(old, new, original, modified)
    clusts_old = [CartesianIndex{2}[]]
    clusts_new=  [CartesianIndex{2}[]]
    for clust in old; push!(clusts_old, findall(x->x==clust, original)); end
    for clust in new; push!(clusts_new, findall(x->x==clust, modified)); end
    popfirst!(clusts_old)
    popfirst!(clusts_new)
    #println(clusts_old, "\n", clusts_new)
    return clusts_old, clusts_new
end


"""
Compute the update to S2 count. Find the contribution from
all white pixel's distance from wpick and bpick, and the difference between them
"""
function update_S2_BN(BN, pix, wpick, bpick, philen, maxrng)
    wp = pix[1][wpick]
    bp = pix[2][bpick]
    F = Tuple.(pix[1])
    x, y= Float32.(first.(F)), Float32.(last.(F))
    dist1 = zeros(Float32, philen) #this can depend on philen
    dist2 = zeros(Float32, philen)

    #pix[2][bpick] = pix[1][wpick]# we need to actually modify x and y, not pix
    x[wpick],y[wpick] = bp[1], bp[2]
    add = fit(Histogram, fragment!(dist1, x, y, wp[1], wp[2]), 0:1.0:maxrng).weights
    x[wpick],y[wpick] = wp[1], wp[2]
    #pix[1][wpick] = pix[2][bpick]
    subtract = fit(Histogram, fragment!(dist2, x, y, bp[1], bp[2]), 0:1.0:maxrng).weights
    #println("ash is of length ", length(subtract))
    return (add .- subtract) .* 2
end

"""
Find the S2 count for an image
"""
function S2_initialise(pix, maxrng)
     (2 .* (blas_stat5(pix[1], maxrng)))
end

"""
Normalise S2 count into S2 proper by dividing by BN and ϕ
"""
function S2_finalise(SN, BN, im)
    dims = size(im)
    BN =  (BN ./ sum(im))
    S2 = ((sum(im)/length(im)) .* (BN ./ (SN ./ length(im))))
end


"""
Loads and thresholds an image
"""
function loadim(fname::String)::Array{Bool, 2}
    img::Array{Bool, 2} = (Gray.(load(fname))) .< 0.50
end
