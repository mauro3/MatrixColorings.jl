using MatrixColorings
using Ragged
using Base.Test

# g00: nothing
adj = RaggedArray(Vector{Int}[])
g00 = Graph(1:0, adj)
# g0: 1
adj = RaggedArray(Vector{Int}[[]])
g0 = Graph(1:1, adj)

# g01: 1 -- 2
adj = RaggedArray(Vector{Int}[[2],[1]])
g01 = Graph(1:2, adj)

# g1:
# 1 -- 2 -- 3 -- 5
#    \ | /
#      4
adj = RaggedArray(Vector{Int}[[2,4], [1,3,4], [2,4,5], [1,2,3], [3]])
g1 = Graph(1:5, adj)

# g2:
# 1 -- 2    3 -- 5
#    \ | 
#      4
adj = RaggedArray(Vector{Int}[[2,4], [1,4], [5], [1,2], [3]])
g2 = Graph(1:5, adj)

# gn
n = 1000
adj = Vector{Int}[]
for i=1:n
    len = rand(1:n-1)
    push!(adj, randperm(n)[1:len])
end
adj = RaggedArray(adj)
gn = Graph(1:n, adj)


# test iteration over nearest neighbors
g,i = 1,1
for g in [g00, g0, g01, g1, g2]#, gn] #uncommenting leads to an error!?
    for i =1:length(g)
        @test Set(g.adjlist[:,i]) == Set(collect(Neigh1(g, i)))
    end
end

@test collect(Neigh2(g0, 1))==Any[]
@test collect(Neigh2(g01, 1))==Any[]

# # test iteration over 2nd neighbors
for g in [g1, g2]
    for v=1:length(g)
        nset = Int[]
        for n in Neigh1(g, v)
            append!(nset, g.adjlist[:,n])
        end
        nset = Set(nset)
        setdiff!(nset, v)
        @test Set{Int}(collect(Neigh2(g, v))) == nset
    end
end

###
# Jacobians
###

function check_color_jac(color, A, dim)
    # Checks whether a coloring col is valid for A along dimension dim
    # Takes a bit of time for large A.
    dim2 = dim==1 ? 2 : 1
    out = zeros(eltype(A), maximum(color), size(A,dim2))
    if dim2==1
        for (j,c) in enumerate(color)
            for jj=1:size(A,1)
                out[c,jj] +=  A[jj,j]
            end
        end
    else
        for (j,c) in enumerate(color)
            for jj=1:size(A,2)
                out[c,jj] +=  A[j,jj]
            end
        end
    end
    @test maximum(out)<=1
    out
end
if VERSION < v"0.4.0-dev"
    # TODO: remove these two for 0.4
    rowvals(S::SparseMatrixCSC) = S.rowval
    nzrange(S::SparseMatrixCSC, col::Integer) = S.colptr[col]:(S.colptr[col+1]-1)
end
function check_color_jac(color, A::SparseMatrixCSC, dim)
    # Checks whether a coloring col is valid for A along dimension dim
    # Takes a bit of time for large A.
    if dim==1
        A = A.'
    end
    out = zeros(eltype(A), maximum(color), size(A,1))
    for (j,c) in enumerate(color)
        for jj in nzrange(A, j)
                out[c,rowvals(A)[jj]] += nonzeros(A)[jj]
        end
    end
    @test maximum(out)<=1
    out
end

# generate graph:
A = [0.0374053  0.0       0.1
     0.0        0.140752  0.380924]
adj = RaggedArray(Vector{Int}[[1+2, 3+2], [2+2,3+2],
                              [1], [2], [1,2]])
cols = Vector{Int}[[1,2], [1,1,2]]
bi = MatrixColorings.matrix2bipartie_graph(A)
@test bi.v1[end] == size(A,1)
@test length(bi.v1)==size(A,1)
@test length(bi.v2)==size(A,2)
@test bi.adjlist==adj
# color graph
col,i = 1,1
for dim=1:2
    col = MatrixColorings.partial_dist2_coloring(bi,dim)
    @test length(col)==size(A,dim)
    @test col==cols[dim]
    check_color_jac(col, A, dim)
end

# generate graph:
A = [0.0374053  1.0       0.1
     1.0        0.140752  0.380924]
adj = RaggedArray(Vector{Int}[[1+2, 2+2, 3+2], [1+2,2+2,3+2],
                              [1,2], [1,2], [1,2]])
cols = Vector{Int}[[1,2], [1,2,3]]
bi = MatrixColorings.matrix2bipartie_graph(A)
@test bi.v1[end] == size(A,1)
@test length(bi.v1)==size(A,1)
@test length(bi.v2)==size(A,2)
@test bi.adjlist==adj
# color graph
col,i = 1,1
for dim=1:2
    col = MatrixColorings.partial_dist2_coloring(bi,dim)
    @test length(col)==size(A,dim)
    @test col==cols[dim]
    check_color_jac(col, A, dim)
end

# generate graph
bi = MatrixColorings.matrix2bipartie_graph(sprand(1,1+1, 0.001)) # warmup

n = 10^4
A = sprand(n,n+1, 0.001)
A.nzval[A.nzval.>0] = 1.0
println("warming up @time:")
@time 1
println("@time for constructing a bipartite graph from $(n)x$n+1 sparse matrix:")
@time bi = MatrixColorings.matrix2bipartie_graph(A)
@test bi.v1[end] == size(A,1)
@test length(bi.v1)==size(A,1)
@test length(bi.v2)==size(A,2)
# color graph
out = 1
j,jj,c,dim = 1,1,1,1
for dim=1:2
    dim2 = dim==1 ? 2 : 1
    println("@time for coloring that graph $(["row","columns"][dim])-wise:")
    @time col = MatrixColorings.partial_dist2_coloring(bi,dim)
    @test length(col)==size(A,dim)
    # sum all columns of same color, no entry should be >1
    check_color_jac(col, A, dim)
end

