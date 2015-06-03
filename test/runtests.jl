using MatrixColorings
using Ragged
using Base.Test

# g00: nothing
adj = RaggedArray(Vector{Int}[])
g00 = Graph{:any}(1:0, adj)
# g0: 1
adj = RaggedArray(Vector{Int}[[]])
g0 = Graph{:any}(1:1, adj)

# g01: 1 -- 2
adj = RaggedArray(Vector{Int}[[2],[1]])
g01 = Graph{:any}(1:2, adj)

# g1:
# 1 -- 2 -- 3 -- 5
#    \ | /
#      4
adj = RaggedArray(Vector{Int}[[2,4], [1,3,4], [2,4,5], [1,2,3], [3]])
g1 = Graph{:any}(1:5, adj)

# g2:
# 1 -- 2    3 -- 5
#    \ | 
#      4
adj = RaggedArray(Vector{Int}[[2,4], [1,4], [5], [1,2], [3]])
g2 = Graph{:any}(1:5, adj)

# gn
n = 1000
adj = Vector{Int}[]
for i=1:n
    len = rand(1:n-1)
    push!(adj, randperm(n)[1:len])
end
adj = RaggedArray(adj)
gn = Graph{:any}(1:n, adj)


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

