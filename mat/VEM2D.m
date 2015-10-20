function [U, hG] = VEM2D(G, f, bc)
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

[S, b, hG] = VEM2D_glob(G, f);

[bcDof, bBC] = VEM2D_bc(G, bc);
b(bcDof == 1) = bBC(bcDof == 1);
SBC = spdiags(ones(Ndof,1),0,Ndof,Ndof);
S(bcDof == 1,:) = SBC(bcDof == 1,:);
b(bcDof == 2) = b(bcDof == 2) + bBC(bcDof == 2);

U = S\b;

end