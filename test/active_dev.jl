using ReconTPCF, Plots


dims, C2, S2, philen = get_C2_S2("test/uniform_circles3.png")
guess, S2n, C2n, S2_BN1, C2_BN1, SN1 = histrecon(dims, C2, S2, philen);
#heatmap(guess, color=:grays, aspect_ratio=1)
save("test38.png", colorview(Gray, guess))
#imshow(guess)
