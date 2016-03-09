function [U, hG] = VEM2D(G, f, bc)
%--------------------------------------------------------------------------
%   Solves the 2D Poisson equation
%
%   -\Delta u = f, x \in \Omega
%           u = g, x \in \partial \Omega
%
%   using VEM.
%
%   G:  MRST Grid
%   f:  Source term
%   bc:  Struct of boundary conditions. See VEM2D_bc for usage.
%--------------------------------------------------------------------------

Nc = G.cells.num;
Ne = G.faces.num;
Nn = G.nodes.num;
Ndof = Nn + Ne + Nc;

[S, b, hG] = VEM2D_glob(G, f);

k = 2;

[bcDof, bBC] = VEM2D_bc(G, bc, k);
b(bcDof == 1) = bBC(bcDof == 1);
SBC = spdiags(ones(Ndof,1),0,Ndof,Ndof);
S(bcDof == 1,:) = SBC(bcDof == 1,:);
b(bcDof == 2) = b(bcDof == 2) + bBC(bcDof == 2);

U = S\b;

end