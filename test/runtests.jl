using ReconTPCF
using Test

ReconTPCF.loadim("C:/Users/James/.julia/dev/ReconTPCF/test/uniform_circles.png")
path = "C:/Users/James/.julia/dev/ReconTPCF"

begin
    test = loadim("uniform_circles.png")
    dims = size(test)#(200,100)
    maxrng = Int64(round(sqrt(dims[1]^2 + dims[2]^2)+1))
    philen = sum(test)#Int(round(0.3*dims[1]*dims[2]))

    SN = SN_comp(dims, maxrng) #Is SN_comp supposed to be multiplied by 2?


    idx = CartesianIndices(dims)
    wpix = findall(x->x==true, test)
    bpix = setdiff(idx, wpix)
    pix = (wpix, bpix)
    S2N = S2_initialise(pix, maxrng)
    S2 = naninf(S2_finalise(SN, S2N, test))
    #S2[length(S2)÷2:end] .= 0
    C2N = C2_initialise(test, maxrng)
    C2 = naninf(S2_finalise(SN, C2N, test))
end

@testset "ReconTPCF.jl" begin
    # Write your tests here.
#ReconTPCF.loadim("uniform_circles.png")
end
