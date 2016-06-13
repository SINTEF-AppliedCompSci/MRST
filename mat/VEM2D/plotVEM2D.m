function h = plotVEM2D(G, sol, k, varargin)
%   Plots the 2D solution obtained from VEM.
%
%   SYNOPSIS:
%       h = plotVEM2D(G, sol, k, varargin)
%
%   REQUIRED PARAMETERS:
%       G       - 2D MRST grid, with sorted edges, G = sortEdges(G),
%                 computed VEM geometry, G = computeVEMGeometry(G).
%       sol     - Solution struct obtained from kth order VEM.
%       k       - Method order. Can be lower than order of VEM used to
%                 obtain sol, but will then use only nodeValues.
%
%   OPTIONAL PARAMETERS:
%       Supports the same optimal parameters as the MATLAB function fill3.
%
%   RETURNS:
%       h       - Figure handle.
%-----------------------------------------------------------------ØSK-2016-

%{
   Copyright (C) 2016 Øystein Strengehagen Klemetsdal. See COPYRIGHT.txt
   for details.
%}

nK = G.cells.num;

for K = 1:nK
    nodeNum = G.cells.nodePos(K):G.cells.nodePos(K+1)-1;
    nodes = G.cells.nodes(nodeNum);
    if k == 1    
        X = G.nodes.coords(nodes,:);
        Z = sol.nodeValues(nodes);
        fill3(X(:,1), X(:,2), Z, Z, varargin{:});
    elseif k == 2
        edgeNum = G.cells.facePos(K):G.cells.facePos(K+1)-1;
        edges = G.cells.faces(edgeNum);
        nN = numel(nodes); nE = numel(edges);
        X = zeros(nN + nE,2);
        X(1:2:end,:) = G.nodes.coords(nodes,:);
        X(2:2:end,:) = G.faces.centroids(edges,:);
        Z = zeros(nN + nE,1);
        Z(1:2:end) = sol.nodeValues(nodes);
        Z(2:2:end) = sol.edgeValues(edges);
        fill3(X(:,1), X(:,2), Z, Z, varargin{:});
    end
    hold on;
    h = gcf;
end
    
        