

include("CorrFuncs.jl")



function histrecon_u(dims, C2, S2, philen)
    guess, pix = make_rand_im(philen, dims)
    display(heatmap(guess, color=:grays, aspect_ratio=1))
    maxrng = Int64(round(sqrt(dims[1]^2 + dims[2]^2)+1))
    wpick = rand(1:length(pix[1]))
    bpick = rand(1:length(pix[2]))
    kern = CartesianIndices((-1:1, -1:1))

    original = label_components(guess, trues(3,3))
    SN = SN_comp(dims, maxrng)
    typeof(SN)
    S2_BN = S2_initialise(pix, maxrng)
    C2_BN = C2_initialise(guess, maxrng)
    S2n = naninf(S2_finalise(SN, S2_BN, guess))
    C2n = naninf(S2_finalise(SN, C2_BN, guess))
    #print(SN, "\n\n", S2_BN, "\n\n", C2_BN)
    Err1 = 2.0
    count = 1
    cont_break=0;
    T0 = 1.0
    cont_surf_opt=0
    pflag = 0
    true_count = 0
    while Err1 > 0.0001
        if (cont_surf_opt !== 1)
            wpick, bpick = pick_rand_coord(pix)
            #println("stoachastic mode")
        else
            surf_res = surf_opt(guess)
            shortlist = surf_rand(surf_res, pix)
            wpick, bpick = pick_rand_coord(shortlist)
            wpick, bpick = find_equivalent(shortlist, pix, (wpick, bpick))
            #println("surface optimisation mode")
            #println("recon shortlist: \twlength = ", length(shortlist[1]), " blength = ", length(shortlist[2]))
            #println("wpick = ", wpick, "\tbpick = ", bpick)


        end

        pix  = swap_pix(guess, pix, wpick, bpick)
        modified = label_components(guess, trues(3,3))
        S2_BN_d = update_S2_BN(S2_BN, pix, wpick, bpick, philen, maxrng)
        C2_BN_d = update_C2_BN(C2_BN, pix, original, modified, wpick, bpick, maxrng)
        #println(length(S2_BN), "\n", length(S2_BN_d))
        S2_BN .+= S2_BN_d
        C2_BN .+= C2_BN_d
        #println("\n", S2_BN_d, "\n")
        #println("C2_BN_d = ", C2_BN_d)
        S2n = naninf(S2_finalise(SN, S2_BN, guess))
        C2n = naninf(S2_finalise(SN, C2_BN, guess))
        #println("C2n = ", C2n)

        true_count +=1

        Err2 =  ediff_C2(C2n, C2) + (ediff_S2(S2n, S2))# +
        #SN_comp(size(guess), length(guess), maxrng)
        delErr = Err2 - Err1
        if count == 1
            T0=-abs(delErr)/log(0.5)
            println("T0 = ", T0)
        end
        T = T0*0.99^count
        if (rand(0:1) < min(exp(-delErr/T), 1))
            count += 1
            cont_break=0
            #print("change commited at count ", count, "\t")
            original = modified
            Err1 = Err2
        else
            swap_pix(guess, pix, wpick, bpick)
            S2_BN .-= S2_BN_d
            C2_BN .-= C2_BN_d
            cont_break += 1

        end

         # we bake in the change here - this must be undone
        if (mod(true_count, 1000) == 0) & (count !=pflag)
            println("\ncount is ", count, "\t  True count is ", true_count, "\t prob of rejection =", min(exp(delErr/T), 1), "\t  S2 Err = ", ediff_S2_fft(S2n,  S2)*1250/length(guess), "\t C2 Err = ",  ediff_C2(C2n, C2))
            println("S2_delerr = ", ediff_S2(S2n, S2), "\tC2_delerr = ", ediff_C2(C2n, C2))
            println("C2[1] = ", C2n[1], "\t\tS2[1] = ", S2n[1])
            if cont_surf_opt==1
                println("Surface optimisation mode")
            else
                println("Stoachastic mode")
            end
            #println(Err1)
            pflag = count
        end



        if count==4000
            cont_surf_opt=1
            #println("Changed to surface optimisation mode")
        end

        if true_count >= 10000

            break
        end




    end
    return guess, S2n, C2n, S2_BN, C2_BN, SN
end




begin
    test = loadim("uniform_circles3.png")
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

guess, S2n, C2n, S2_BN1, C2_BN1, SN1 = histrecon_u(dims, C2, S2, philen);
heatmap(guess, color=:grays, aspect_ratio=1)
save("test38.png", colorview(Gray, guess))
imshow(guess)




#dims = (500,500)
#philen = Int(round(0.3*dims[1]*dims[1]))

#@time SN, BN = C2_init_partial(f, 100)
#@time subtract, add = update_C2_BN(BN, pix, original, modified, wpick, bpick)
#C2n, C2_BN = C2_pdist4(f, 100, SN)
a = [0]
dims = (5,5)
maxrng = Int64(round(sqrt(dims[1]^2 + dims[2]^2)+1))
philen = Int(round(0.3*dims[1]*dims[2]))
SN = SN_comp(dims, maxrng)

@time blas_stat5(vec(CartesianIndices(dims)), maxrng)


f, pix = make_rand_im(philen, dims)
heatmap(f, color=:grays, aspect_ratio=1)

    #pix = (wpix, bpix)
wpick = rand(1:length(pix[1]))
bpick = rand(1:length(pix[2]))
kern = CartesianIndices((-1:1, -1:1))

original = label_components(f, trues(3,3))
pix = swap_pix(f, pix, wpick, bpick)
modified = label_components(f, trues(3,3))
#a = pix


S2_BN = S2_initialise(pix, maxrng)
C2_BN = C2_initialise(f, maxrng)
S2n3 = S2_finalise(SN, S2_BN, f)
S2_finalise(SN, C2_BN, f)
S2_BN_d = update_S2_BN(S2_BN, pix, wpick, bpick, philen, maxrng)
C2_BN_d = update_C2_BN(C2_BN, pix, original, modified, wpick, bpick, maxrng)
push!(a, C2_BN_d[1])
blas_stat4(pix[1], 1.0, 100)

original = label_components(f, trues(3,3))
pix = swap_pix(f, pix, wpick, bpick)
modified = label_components(f, trues(3,3))
q = add - subtract
q = add - subtract
BN .= BN .+ q
#c1 = copy(change)
q = [0.0, 2.236068, 1.4142135, 1.0, 2.828427, 2.236068, 2.0, 0.0, 0.0, 1.0, 1.4142135, 1.0, 1.4142135, 0.0, 1.0, 2.0, 0.0, 2.236068, 1.4142135, 1.0, 0.0, 1.0, 1.4142135, 1.0, 1.4142135, 2.236068, 2.0, 2.236068, 0.0, 1.0, 0.0, 1.0, 2.0, 1.0, 1.4142135, 2.236068, 0.0, 1.0, 2.0, 1.0, 1.4142135, 2.236068, 2.0, 2.236068, 2.828427]
ash(q, rng=0:1.0:(maxrng-1), m =1).counts

@time C2_pdist4(f, 100, SN)
heatmap(original, color=:grays, aspect_ratio=1)
heatmap(modified, color=:grays, aspect_ratio=1)
imshow(original)
imshow(modified)
CartesianIndex(1,1)

a = collect(0:4).^2
b = a'
fit(Histogram, vec(sqrt.(a.+b)), 0:1.0:(maxrng)).weights
