function sol = cellAverages(G, sol)

nK = G.cells.num;

nNK = diff(G.cells.nodePos);

cellAverages = zeros(nK, 1);
for K = 1:nK
    nodeNum = G.cells.nodePos(K):G.cells.nodePos(K+1)-1;
    nodes = G.cells.nodes(nodeNum);
    cellAverages(K) = sum(sol.nodeValues(nodes))/nNK(K);
end
sol.('cellAverages') = cellAverages;
    