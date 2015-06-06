module MatrixColorings
if VERSION < v"0.4.0-dev"
    using Docile
end

using Ragged

import Base: start, next, done, length

export Graph, Neigh1, Neigh2

# this follows "ColPack: Software for Graph Coloring and Related
# Problems in Scientific Computing" GEBREMEDHIN & et al

typealias MatrixLike Any # maybe have a trait here one day

## Misc helper functions
#####

# Allocates all vectors inside a Vector{Vector{T}} and element type
# eltype(eltype(vec)).
allocate!{T}(vec::Vector{T}, n::Int) =
    (vec[:] = Vector{T}[ Array(eltype(T), n) for s=1:length(vec)])

## Graphs
#########
# Reinventing the wheel?

# Compact storage for an undirected graph.
# Neighbors of vertex i are adjlist[i]
abstract AGraph
# assumes fields
# v::UnitRange{Int}  # vertices
# adjlist::RaggedArray{Int,Int} # adjacency list, i.e. adjlist[v]
                                # gives connected vertices.
length(g::AGraph) = length(g.v)

immutable Graph <: AGraph
    v::UnitRange{Int}  # vertices
    adjlist::RaggedArray{Int,Int} # adjacency list
    work::BitArray{1} # work array for various algos
    function Graph(v,a)
        @assert size(a,2)==length(v)
        length(v)>0 && @assert v[1]==1
        new(v,a,falses(length(v)))
    end
end
immutable BiGraph <: AGraph
    v::UnitRange{Int}  # vertices
    adjlist::RaggedArray{Int,Int} # adjacency list
    work::BitArray{1} # work array for various algos
    v1::UnitRange{Int} # first partition:
    v2::UnitRange{Int} # second partition
    function BiGraph(v,a,v1,v2)
        @assert size(a,2)==length(v)
        if length(v)>0
            @assert v[1]==1
            @assert v[1]==1
            @assert v1[1]==1
            @assert v1[end]+1==v2[1]
            @assert v2[end]==v[end]
        else
            @assert length(v1)==length(v2)==0
        end
        new(v,a,falses(length(v)),v1,v2)
    end
end


@doc "Iterates over all direct neighbor vertices." -> 
immutable Neigh1{G<:AGraph}
    g::G
    v::Int
end
start(n1::Neigh1) = 1 # state is row-index of adjlist
next(n1::Neigh1, state) = (n1.g.adjlist[state,n1.v], state+1)
done(n1::Neigh1, state) = Ragged.collength(n1.g.adjlist,n1.v)+1==state

@doc "Iterates over all 2nd neighbor vertices." -> 
immutable Neigh2{G<:AGraph}
    g::G
    v::Int
end
# The problem here is that we don't know we've reached the end until
# we've reached it, thus use the dict.jl skip_deleted trick.

# Type to hold the state:
type Neigh2State{G<:AGraph}
    n1::Neigh1{G}  # the current Neigh1 to iterate over
    n1_state::Int  # the state of n1
    n2_state::Int  # the state encoding which n1 we're on
    item::Int      # the actual item
end
# Function which does all the logic
function _get_next(n2::Neigh2, s::Neigh2State)
    # Iterate over n1, discarding an item if it has been returned
    # already or if it is n2.v.
    visited = n2.g.work
    num_n1 = Ragged.collength(n2.g.adjlist, n2.v)
    @inbounds while visited[s.item] # if item has been visited, try next
        while done(s.n1, s.n1_state) # if n1 is done go onto next first neighbor
            s.n2_state += 1 # onto the next neigh
            if s.n2_state <= num_n1
                s.n1 = Neigh1(n2.g, n2.g.adjlist[s.n2_state, n2.v])
                s.n1_state = start(s.n1)
            else
                s.n2_state = -1  # this indicates completely done,
                                 # i.e. no next item.
                return s
            end
        end
        s.item, s.n1_state = next(s.n1, s.n1_state)
    end
    visited[s.item] = true
    return s
end    

# State is (Neigh1, state-of-Neigh1, i in g.adjlist[i,v])
function start(n2::Neigh2)
    if length(n2.g)==0 || Ragged.collength(n2.g.adjlist,n2.v)==0
        # no 1st neighbors
        n1 = Neigh1(n2.g, -1)
        n1_state = -1
        n2_state = -1
        return Neigh2State(n1, n1_state, n2_state, 0)
    else
        fill!(n2.g.work, false) # initialize work array to keep track
                                # of visited vertices.
        n2.g.work[n2.v] = true  # 
        item = n2.v
        n1 = Neigh1(n2.g, n2.g.adjlist[1,n2.v])
        n1_state, n2_state = start(n1), 1
        return _get_next(n2, Neigh2State(n1, n1_state, n2_state, item))
    end
end
next(n2::Neigh2, state) = (state.item, _get_next(n2, state))
done(n2::Neigh2, state) = state.n2_state==-1

# function Neigh2(g::AGraph, v)
#     # stopgap until above is fixed
#     nset = Int[]
#     for n in Neigh1(g, v)
#         append!(nset, g.adjlist[:,n])
#     end
#     nset = Set(nset)
#     setdiff!(nset, v)
#     return nset
# end

# For Jacobians use "partial distance-2 coloring"

@doc """

\"The structure of a Jacobian matrix A is represented by the bipartite
graph Gb(A) = (V₁, V₂, E), where the vertex sets V₁ and V₂ represent
the rows and columns of A, respectively, and each nonzero matrix entry
A_ ij is represented by the edge (r_i , c_j) in E.\"

""" ->
function gen_bipartie_graph(A::MatrixLike)
    nr,nc = size(A)
    nv = nr+nc
    part = nr+1
    v = 1:nv
    v1 = 1:nr
    v2 = nr+1:nv
    adj = Array(Vector{Int}, nv)
    allocate!(adj, 0)
    for j=1:nc,i=1:nr
        if A[i,j]!=0
            vv1, vv2 = v1[i], v2[j]
            push!(adj[vv1], vv2)
            push!(adj[vv2], vv1)
        end
    end
    return BiGraph(v, RaggedArray(adj), v1, v2)
end

@doc """I think sorting can improve the coloring, but no idea how to
sort."""->
function sort_vertices!(g::AGraph)
    nothing
    # probably need to return a permeation to reverse the sort.
end

""" 
\"In a bipartite graph G b = (V 1 , V 2 , E), a partial distance-2
coloring on the vertex set V 1 (or V 2 ) is an assignment of colors to
the vertices in V 1 (or V 2 ) such that a pair of vertices connected
by a path of length exactly two edges receives different colors. The
term partial, which is sometimes omitted when the context is clear, is
used to emphasize that the other vertex set remains uncolored.\" 
"""
function partial_dist2_coloring(g::BiGraph, part=1)
    vv = part==1 ? g.v1 : g.v2
    #sort_vertices!(g::AGraph)
    nvv = length(vv)

    # This is the output.  Possible colors range from 1:nnv, indexed
    # by index into vv.
    color = zeros(Int, nvv)
    # A work array holding vertex numbers, indexed by
    # color. Initialize with a non-vertex number:
    forbiddenColors = zeros(Int,nvv)
    # - forbiddenColors[c] == v indicates that the color c is
    #   impermissible for the vertex v
    # - forbiddenColors[c] != v indicates that c is a candidate color
    #   for the vertex v, regardless of the actual value in
    #   forbiddenColors[c].
    for (i,v) in enumerate(vv)
        for w in Neigh2(g, v) # replace here for other types of coloring
            ii = w-vv[1]+1 # index of w into vv
            @assert w==vv[ii] # TODO remove
            if color[ii]>0
                forbiddenColors[color[ii]] = i
            end
        end
        for (col,fc) in enumerate(forbiddenColors)
            if fc!=i
                color[i] = col
                break
            end
        end
    end    
    return color
end

# 4.2. Algorithms for coloring problems on bipartite graphs

end # module
