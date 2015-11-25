function p = plotGridWithDofs(G,bc)

Nc = G.cells.num;
Ne = G.faces.num;
Nn = G.nodes.num;
Ndof = Nn + Ne + Nc;

% boundaryNodes = zeros(Ndof,1);
% neighbors = G.faces.neighbors;
% for e = 1:Ne
%     if neighbors(e,1) == 0 || neighbors(e,2) == 0
%         nodeNum = G.faces.nodePos(e):G.faces.nodePos(e+1)-1;
%         boundaryNodes([G.faces.nodes(nodeNum)', Nn + e]) = 1;
%     end
% end
% 
% boundaryNodes = find(boundaryNodes);

X = [G.nodes.coords ; G.faces.centroids ; G.cells.centroids];

plotGrid(G)

hold on
bcFaces = bc.bcFaces;
Nbc = numel(bcFaces);
for b = 1:Nbc
    bcEdges = bcFaces{b};
    Ne = numel(bcEdges);
    bcNodes = zeros(1,2*Ne);
    for e = 1:Ne
        nodeNum = G.faces.nodePos(bcEdges(e)):G.faces.nodePos(bcEdges(e)+1)-1;
        nodes = G.faces.nodes(nodeNum);
        nodes = fixDim(nodes);
        bcNodes(3*e-2:3*e) = [nodes, Nn + bcEdges(e)];
    end
    plot(X(bcNodes,1),X(bcNodes,2),'o');
end  
plot(X(:,1),X(:,2),'.k');
hold off
axis([0 1 0 1]);
xlabel('x'); ylabel('y');
legend('Cells', 'Boundary dofs', 'Internal dofs');

end