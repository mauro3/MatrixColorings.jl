# MatrixColorings

This functionality is now included in the DiffEq-universe, in particular https://github.com/JuliaDiffEq/SparseDiffTools.jl.

[![Build Status](https://travis-ci.org/mauro3/MatrixColorings.jl.svg?branch=master)](https://travis-ci.org/mauro3/MatrixColorings.jl)

This implements coloring algorithms for matrices.  Should be useful
for finite difference or automatic differentiation calculations of
Jacobians and Hessians.

Uses Ragged.jl, which is on BitBucket, add with
```
Pkg.clone("git@bitbucket.org:maurow/ragged.jl.git", "Ragged")
```

Example of numerical Jacobian (see `example/jacobian_eval.jl`):

```julia
using MatrixColorings
jac = ... # sparse Jacobian allocation
S = color_jac(jac) # seed matrix based on sparsity structure of jac
                   # such that dense_jac = jac * S

function numjac_c!(t, y, jac, S)
    # Finite difference Jacobian, updates jac in-place.
    nc = size(S,2) # number of colors
    
    edelta  = eps(4.0)*1e7 # the delta for the finite differences
    dof = length(y)
    f0 = similar(y)   # f(d)
    fn!(t, y, f0)
    ddy = deepcopy(y) # d+dy
    f1 = similar(y)   # f(d+dy)
    for c in 1:nc # loop over all colors
        add_delta!(ddy, S, c, edelta)  # adds a delta to all indices of color c
        fn!(t, ddy, f1) # f(d+dy)
        for j=1:dof
            f1[j] = (f1[j]-f0[j])/edelta
        end
        recover_jac!(jac, f1, S, c) # recovers the columns of jac
                                    # which correspond to color c
        remove_delta!(ddy, S, c, edelta) # remove delta from ddy
    end
    return nothing
end
numjac!(0, y, jac, S) # calculates Jacobian
```

For the example Jacobian of 20,000x20,000 `example/jacobian_eval.jl`,
this method is 3000x faster than the brute force finite difference and
only 10x slower than the analytic one.

# Other packages
- [ReverseDiffSparse.jl](https://github.com/mlubin/ReverseDiffSparse.jl/blob/master/src/coloring.jl)
  has coloring algorithms for Hessians using acyclic coloring.

# References
- [Hossain, S. and Steihaug, T. 2012](http://www.tandfonline.com/doi/abs/10.1080/10556788.2012.693927)
  might be state of the art for Jacobians
