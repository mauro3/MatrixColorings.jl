# Example evaluating a numerical a Jacobian.  Uses the Brusselator in IVPtestsuite.jl

# using PyPlot
using MatrixColorings

const N = 10^4 # active gridpoints
dof = 2N # degrees of freedom (boundary points are not free)
dx = 1/(N+1)
x = dx:dx:1-dx
T = Float64

const alpha = 1/50
const gamma = alpha/dx^2

# BC
const one_ = one(T)
const three = 3*one_
ubc1(t) = one_
ubc2(t) = one_
vbc1(t) = three
vbc2(t) = three

# immutable Pair
#     i::Int
#     j::Int
# end
# inds(i) = Pair(2i-1, 2i)
inds(i) = (2i-1, 2i)
function fn!(t,y, dydt)
    # the ode function dydt = f(t,y)
    #
    # Note the u and v are staggered:
    # u = y[1:2:end], v=y[2:2:end]
    
    # setup for first step
    i = 1
    iu, iv = inds(i)

    ui0 = ubc1(t); vi0 = vbc1(t) # BC
    ui1 = y[iu];   vi1 = y[iv]
    ui2 = y[iu+2]; vi2 = y[iv+2]
    while i<=N
        dydt[iu] = 1 + ui1^2*vi1 - 4*ui1 + gamma*(ui0 - 2*ui1 + ui2)
        dydt[iv] = 3*ui1 - ui1^2*vi1     + gamma*(vi0 - 2*vi1 + vi2)
        # setup for next step:
        i += 1
        iu, iv = inds(i)
        ui0 = ui1; vi0 = vi1
        ui1 = ui2; vi1 = vi2
        if i<N
            ui2 = y[iu+2]
            vi2 = y[iv+2]
        else # on last step do BC
            ui2 = ubc2(t)
            vi2 = vbc2(t)
        end
    end
    return nothing
end
fn!(;T_::Type=T, dof_=dof) = zeros(T_,dof_) # returns array for storage for dydt

# Analytic Jacobian: it is a banded matrix with upper and lower
# bandwidth 2, c.f. W&H p.148, but Julia does not support banded
# matrices yet.  Thus use a sparse matrix.
function jac!(t,y,dfdy)
    # The Jacobian of fn
    ii = 1 # direct index into dfdy.nzval
    for i=1:N # iterate over double-columns
        iu, iv = inds(i)
        ui = y[iu]
        vi = y[iv]
        # iu: column belonging to d/du(i)
        if iu-2>=1
            # du'(i-1)/du(i)
            dfdy.nzval[ii] = gamma
            ii +=1
        end
        if iv-2>=1
            # dv'(i-1)/du(i)==0
            dfdy.nzval[ii] = 0
            ii +=1
        end
        # du'(i)/du(i)
        dfdy.nzval[ii]   = 2*ui*vi - 4 -2*gamma
        ii += 1
        # dv'(i)/du(i)
        dfdy.nzval[ii] = 3 - 2*ui*vi
        ii +=1
        if  iu+2<=2N
            # du'(i+1)/du(i)
            dfdy.nzval[ii] = gamma
            ii +=1
        end
        
        # iv: column belonging to d/dv(i)
        if iv-2>=1
            # dv'(i-1)/dv(i)
            dfdy.nzval[ii] = gamma
            ii +=1
        end
        # du'(i)/dv(i)
        dfdy.nzval[ii] = ui^2
        ii +=1
        # dv'(i)/dv(i)
        dfdy.nzval[ii]   = -ui^2 -2*gamma
        ii +=1
        if  iu+2<=2N
            # du'(i+1)/dv(i)==0
            dfdy.nzval[ii] = 0
            ii +=1
        end
        if iv+2<=2N
            # dv'(i+1)/dv(i)
            dfdy.nzval[ii] = gamma
            ii +=1
        end
    end
    return nothing
end
function jac!(;T_::Type=T, dof_=dof)
    # returns an initialized sparse matrix
    B = (ones(T_,dof_-2), ones(T_,dof_-1), ones(T_,dof_), ones(T_,dof_-1), ones(T_,dof_-2))        
    spdiagm(B, (-2,-1,0,1,2))
end

####
# Numerical Jacobians
####

# Non-coloring method:
######################
# Copied & adapted from DASSL:
# generate a function that computes approximate jacobian using forward
# finite differences
if VERSION < v"0.4.0-dev"
    # TODO: remove these two for 0.4
    rowvals(S::SparseMatrixCSC) = S.rowval
    nzrange(S::SparseMatrixCSC, col::Integer) = S.colptr[col]:(S.colptr[col+1]-1)
end
function numjac!(t, y, jac)
    # updates jac in-place
        
    edelta  = eps(4.0)*1e7 # trial and error
    dof = length(y)
    f1 = similar(y)
    tmpy = similar(y) # d+dy
    f0 = similar(y)
    fn!(t, y, f0)
    for i=1:length(y)
        copy!(tmpy, y)
        tmpy[i] += edelta
        fn!(t, tmpy, f1)
        for j=1:dof
            tmpy[j] = (f1[j]-f0[j])/edelta
        end
        for nz in nzrange(jac,i)
            nonzeros(jac)[nz] = tmpy[rowvals(jac)[nz]]
        end
        # which, for large matrices, is much faster than:
        # jac[:,i] = (f1-f0)/edelta
    end
    return nothing
end


# Coloring method:
function numjac_c!(t, y, jac, S)
    # updates jac in-place
    nc = size(S,2) # number of colors
    
    edelta  = eps(4.0)*1e7 # trial and error
    dof = length(y)
    f0 = similar(y)
    fn!(t, y, f0)
    tmpy = similar(y) # d+dy
    f1 = similar(y)   # f(d+dy)
    for c in 1:nc
        copy!(tmpy, y)
        for i in nzrange(S,c)
            tmpy[rowvals(S)[i]] += edelta
        end
        fn!(t, tmpy, f1)
        for j=1:dof
            tmpy[j] = (f1[j]-f0[j])/edelta
        end
        recover_jac!(jac, tmpy, S, c)
    end
    return nothing
end

function numjac_c!(t, y, jac)
    # updates jac in-place, uses its sparsity structure for coloring
    S = color_jac(jac) # seed matrix
    numjac_c!(t, y, jac, S)
    return nothing
end

####
# Compare
####

t = 0.0 # not dependent on t
# initialize
dydt = fn!()
dfdy_ana = jac!()
dfdy_num = jac!()
dfdy_num_c = jac!()

y = fn!()
#IC
u0(x) = 1 + 0.5*sin(2π*x)
v0(x) = 3 * ones(length(x))
y[1:2:2N] = u0(x)
y[2:2:2N] = v0(x)
# warm ups:
@time 1
jac!(t, y, dfdy_ana)
numjac!(t, y, dfdy_num)

#ys = (y, rand(T, dof)*3+0.5)
for y in (y,)
    # calculate, time and compare
    fn!(t, y, dydt)
    println("@time for analytic:")
    @time jac!(t, y, dfdy_ana)
    println("@time for numeric with no coloring:")
    @time numjac!(t, y, dfdy_num)
    println("@time for numeric with coloring:")
    @time numjac_c!(t, y, dfdy_num_c)
    S = color_jac(dfdy_num_c) # seed matrix
    println("@time for numeric with coloring with provided seed matrix:")
    @time numjac_c!(t, y, dfdy_num_c,S)


    @show maximum(dfdy_ana-dfdy_num)
    

    @show maximum(dfdy_ana-dfdy_num_c)
end

