function p = plotGridWithDofs(G)

Nc = G.cells.num;
Ne = G.faces.num;
Nn = G.nodes.num;
Ndof = Nn + Ne + Nc;

baricenters = zeros(Nc,2);
for c = 1:Nc
    nodeNum = G.cells.nodePos(c):G.cells.nodePos(c+1)-1;        
    nodes = G.cells.nodes(nodeNum);
    X = G.nodes.coords(nodes,:);
    [~, baricenters(c,:)] = baric(X);
end

boundaryNodes = zeros(Ndof,1);
neighbors = G.faces.neighbors;
for e = 1:Ne
    if neighbors(e,1) == 0 || neighbors(e,2) == 0
        nodeNum = G.faces.nodePos(e):G.faces.nodePos(e+1)-1;
        boundaryNodes([G.faces.nodes(nodeNum)', Nn + e]) = 1;
    end
end

boundaryNodes = find(boundaryNodes);

X = [G.nodes.coords ; G.faces.centroids ; baricenters];

plotGrid(G)
hold on
    plot(X(boundaryNodes,1), X(boundaryNodes,2),'or');
    plot(X(:,1),X(:,2),'.k');
hold off
axis([0 1 0 1]);
xlabel('x'); ylabel('y');
legend('Cells', 'Boundary dofs', 'Internal dofs');

end