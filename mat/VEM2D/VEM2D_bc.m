function [A, b] = VEM2D_bc(G, A, b, bc, k)
%   Incorporates boundary conditions in stiffness matrix A and load term b,
%   obtained using a kth order VEM.
%
%   SYNOPSIS:
%       [A, b] = VEM2D_bc(G, A, b, bc, k)
%
%   DESCRIPTION:
%       Incorporates boundary conditions in stiffness matrix A and load
%       term b, obtained using a kth order VEM. Dirichlet conditions are
%       set directly in the load term, and the corresponding row of A is
%       changed to the corresponding row of the identity matrix. Neumann
%       conditions adds a boundary integral to the load term, which is
%       evaluated using a k+1 accurate quadrature rule. See [1] for
%       details.
%
%   REQUIRED PARAMETERS:
%       G   - MRST grid.
%       A   - Global VEM stiffness matrix.
%       b   - Global VEM load term.
%       bc  - Boundary condition struct constructed using VEM2D_addBC.
%       k   - Method order.
%
%   RETURNS:
%       A   - Global stiffness matrix with boundary condition
%             modifications.
%       b   - Global load term with boundary condition modifications.    
%
%   REFERENCES:
%       [1] - The virtual element method as a common framework for
%             finite element and finite difference methods - Numerical
%             and theoretical analysis.
%-----------------------------------------------------------------ØSK-2016-

%{
   Copyright (C) 2016 Øystein Strengehagen Klemetsdal. See Copyright.txt
   for details.
%}

assert(k == 1 | k == 2); 

NK = G.nodes.num + G.faces.num*(k-1) + G.cells.num*k*(k-1)/2;

edges   = bc.face(strcmp(bc.type,'flux'));
edgeLengths = G.faces.areas(edges);

nodeNum = mcolon(G.faces.nodePos(edges), G.faces.nodePos(edges +1)-1);
nodes   = G.faces.nodes(nodeNum);
nN = numel(nodes);
uNodes = unique(nodes);
nUN = numel(uNodes);
S = (repmat(nodes,1,nUN) == repmat(uNodes',nN,1))';

vals = bc.value(strcmp(bc.type,'flux'),:);

if k == 1
    vals = vals.*[edgeLengths/6, edgeLengths/6, edgeLengths/3];
    vals = bsxfun(@plus, vals(:,1:2), vals(:,3));
    vals = reshape(vals',[],1);
    vals = S*vals;
    dofVec = uNodes';
elseif k == 2
    vals = vals.*[edgeLengths/6, edgeLengths/6, edgeLengths*2/3];
    vals = [S*reshape(vals(:,1:2)',[],1); vals(:,3)];
    dofVec = [uNodes', edges' + G.nodes.num];
end
b(dofVec) = b(dofVec) + vals;

edges   = bc.face(strcmp(bc.type,'pressure'));
nodeNum = mcolon(G.faces.nodePos(edges), G.faces.nodePos(edges +1)-1);
nodes   = G.faces.nodes(nodeNum);
vals = bc.value(strcmp(bc.type,'pressure'),:);

if k == 1
    vals = reshape(vals(:,1:2)',[],1);
    dofVec = nodes';
elseif k == 2
    dofVec = [nodes', edges' + G.nodes.num];
    vals = [reshape(vals(:,1:2)',[],1); vals(:,3)];
end
b(dofVec) = vals;
I = spdiags(ones(NK,1),0,NK,NK);
A(dofVec,:) = I(dofVec,:);

end 