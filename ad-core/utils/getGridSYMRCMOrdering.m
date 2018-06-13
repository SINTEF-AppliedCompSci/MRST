function ordering = getGridSYMRCMOrdering(G)
    N = getNeighbourship(G);
    A = getConnectivityMatrix(N, 'Topological');
    ordering = symrcm(A);
end