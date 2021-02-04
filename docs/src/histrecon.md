# Reconstruction

```@meta
CurrentModule = ReconTPCF
DocTestSetup = quote
    using ReconTPCF
end
```

The reconstruction algorithm uses simulated annealing to drive a system to minimise
a cost function. This cost function is the mismatch in the two-point correlation
function between an ideal state and a present one. After initialising the system
with a trial (random) binary image, we move towards an image with similar ``C_{2}``
and ``S_{2}``. This is done via a series of discrete pixel substitutions. After
each substitution, we calculate how ``C_{2}`` and ``S_{2}`` have changed by subtracting
those pixel's original contributions. We then add the new contribution to ``C_{2}``
and ``S_{2}`` from the substitution. By this we reduce the needed number of
computations.

This is simple for ``S_{2}`` as we simply compute all contributions from the pixels
involved in the swap and find the difference. For ``C_{2}`` we must first determine
which collections of pixels are involved in the interaction, then compute their
contributions to ``C_{2}`` and ``S_{2}``.

## Adjusting the reconstruction algorithm

The temperature is currently set automatically, but this could instead be driven by a
schedule more complex than geometric decay.

Surface optimisation is key to the success of the algorithm as it reduces the set of
pixels permitted to mutate. Despite this, it biases the resultant reconstruction
towards larger particles, and hinders the closing of bridges. Several convolutional
stages could help to reduce this tendency and better settle the system.

It may be worth using more aggressive shuffles to begin with and fully recompute
``C_{2}`` and ``S_{2}``.


```@docs
histrecon(dims, C2, S2, philen)
get_C2_S2(fname)
```
