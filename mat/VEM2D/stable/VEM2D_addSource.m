function q = VEM2D_addSource(G, cells, Q, k)

if k == 1
    
    q = Q./G.cells.