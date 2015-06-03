module MatrixColorings
if VERSION < v"0.4.0-dev"
    using Docile
end

using Ragged

import Base: start, next, done, length

export Graph, Neigh1, Neigh2

# this follows "ColPack: Software for Graph Coloring and Related
# Problems in Scientific Computing" GEBREMEDHIN & et al

## Graphs
#########
# Reinventing the wheel?

# Compact storage for an undirected graph.
# Neighbors of vertex i are adjlist[i]
immutable Graph{Typ}
    v::UnitRange{Int}  # vertices
    adjlist::RaggedArray{Int,Int} # adjacency list, i.e. adjlist[v]
                                  # gives connected vertices.
end
length(g::Graph) = length(g.v)
@doc "Iterates over all direct neighbor vertices." -> 
immutable Neigh1
    g::Graph
    v::Int
end
start(n1::Neigh1) = 1 # state is row-index of adjlist
next(n1::Neigh1, state) = (n1.g.adjlist[state,n1.v], state+1)
done(n1::Neigh1, state) = Ragged.collength(n1.g.adjlist,n1.v)+1==state

@doc "Iterates over all 2nd neighbor vertices." -> 
# immutable Neigh2
#     g::Graph
#     v::Int
# end
# # The problem here is that we don't know we've reached the end until
# # we've reached it, thus use the dict.jl skip_deleted trick.

# # Type to hold the state:
# type Neigh2State
#     n1::Neigh1     # the Neigh1 to iterate over
#     n1_state::Int  # the state of n1
#     n2_state::Int  # the state encoding which n1 we're on
#     visited::BitArray{1} # to keep track of visits
#     item::Int # the actual item
# end
# # Function which does all the logic
# function _get_next(n2::Neigh2, s::Neigh2State)
#     # iterate over n1, discarding item if it has been returned already
#     # or if it is n2.v
#     s.item = n2.v
#     while s.visited[s.item]
#         if done(s.n1, s.n1_state) # if n1 is done go onto next first neighbor
#             s.n2_state += 1 # onto the next neigh
#             if s.n2_state<Ragged.collength(n2.g.adjlist, n2.v)
#                 s.n1 = Neigh1(n2.g, n2.g.adjlist[s.n2_state, n2.v])
#                 s.n1_state = start(s.n1)
#             else
#                 s.n2_state = -1  # this indicates completely done
#                 break
#             end
#         end
#         s.item, s.n1_state = next(s.n1, s.n1_state)
#     end
#     s.visited[s.item] = true
#     return s
# end    

# # State is (Neigh1, state-of-Neigh1, i in g.adjlist[i,v])
# function start(n2::Neigh2)
#     if length(n2.g)==0 || Ragged.collength(n2.g.adjlist,n2.v)==0
#         # no 1st neighbors
#         n1 = Neigh1(n2.g, -1)
#         n1_state = -1
#         n2_state = -1
#         visited = falses(0)
#         return Neigh2State(n1, n1_state, n2_state, visited, 0)
#     else
#         n1 = Neigh1(n2.g, n2.g.adjlist[1,n2.v])
#         n1_state, n2_state = start(n1), 1
#         visited = falses(length(n2.g))
#         visited[n2.v] = true
#         return _get_next(n2, Neigh2State(n1, n1_state, n2_state, visited, 0))
#     end
# end
# next(n2::Neigh2, state) = (state.item, _get_next(n2, state))
# done(n2::Neigh2, state) = state.n2_state==-1
function Neigh2(g::Graph, v)
    # stopgap until above is fixed
    nset = Int[]
    for n in Neigh1(g, v)
        append!(nset, g.adjlist[:,n])
    end
    nset = Set(nset)
    setdiff!(nset, v)
    return nset
end

# For Jacobians use "partial distance-2 coloring"

@doc """

\"The structure of a Jacobian matrix A is represented by the bipartite
graph Gb(A) = (V₁, V₂, E), where the vertex sets V₁ and V₂ represent
the rows and columns of A, respectively, and each nonzero matrix entry
A_ ij is represented by the edge (r_i , c_j) in E.\"

""" ->
function gen_bipartie_graph()

    return out::Graph{:bipartite}
end

@doc """I think sorting can improve the coloring, but no idea how to
sort."""->
function sort_vertices!(g::Graph)
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
function partial_dist2_coloring(g::Graph)
    vv = g.v
    #sort_vertices!(g::Graph)
    nvv = length(vv)

    # This is the output.  Possible colors range from 1:nnv, indexed
    # by vertex number.
    color = zeros(Int, nvv)

    # A work array holding vertex numbers, indexed by
    # color. Initialize with vertex, only need to initialize once.
    forbiddenColors = ones(Int,nvv) + nvv+1
    # - forbiddenColors[c] == v indicates that the color c is
    #   impermissible for the vertex v
    # - forbiddenColors[c] != v indicates that c is a candidate color
    #   for the vertex v, regardless of the actual value in
    #   forbiddenColors[c].
    for v in vv
        for w in neigh2(v) # replace here for other types of coloring
            if color[w]>0
                forbiddenColors[color[w]] = v
            end
        end
        for (i,fc) in enumerate(forbiddenColors)
            if fc!=v
                color[v] = i
                break
            end
        end
    end    
    return color
end

# 4.2. Algorithms for coloring problems on bipartite graphs

end # module
