using MatrixColorings
#using ProfileView

# generate graph
n = 1000
A = sprand(n,n+1, 0.1)
bi = MatrixColorings.matrix2bipartie_graph(A)
bi = MatrixColorings.matrix2bipartie_graph(A)
# color graph
i=1
col = MatrixColorings.partial_dist2_coloring(bi,i)
@time col = MatrixColorings.partial_dist2_coloring(bi,i)
#ProfileView.view()


nothing
