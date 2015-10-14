function U = VEM2D(G, f, g)
%--------------------------------------------------------------------------
%   -\Delta u = f, x \in \Omega
%           u = g, x \in \partial \Omega
%
%   G:  Grid
%   f:  Force term
%   g:  Struct of boundary conditions ...
%--------------------------------------------------------------------------

Nc = G.cells.num;
Ne = G.faces.num;
Nn = G.nodes.num;
Ndof = Nn + Ne + Nc;

[S, b] = VEM2D_glob(G, f);

boundaryNodes = zeros(Ndof,1);
neighbors = G.faces.neighbors;
for e = 1:Ne
    if neighbors(e,1) == 0 || neighbors(e,2) == 0
        nodeNum = G.faces.nodePos(e):G.faces.nodePos(e+1)-1;
        boundaryNodes([G.faces.nodes(nodeNum)', Nn + e]) = 1;
    end
end

boundaryNodes = find(boundaryNodes);

SBC = eye(Ndof);

S(boundaryNodes,:) = SBC(boundaryNodes,:);

baricenters = zeros(Nc,2);
for c = 1:Nc
    nodeNum = G.cells.nodePos(c):G.cells.nodePos(c+1)-1;        
    nodes = G.cells.nodes(nodeNum);
    X = G.nodes.coords(nodes,:);
    [~, baricenters(c,:)] = baric(X);
end

X = [G.nodes.coords ; G.faces.centroids ; baricenters];
b(boundaryNodes) = g(X(boundaryNodes,:));

U = S\b;

end