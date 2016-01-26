function [S,b] = VEM3D_glob(G,f,bc)

nN = G.nodes.num;
nE = G.edges.num;
nF = G.faces.num;
nK = G.cells.num;
k = 2;

N = nN + nE*(k-1) + nF*k*(k-1)/2 + nK*k*(k^2-1)/6;

S = sparse(N,N);
b = sparse(N,1);

for i = 1:nK

    [Sl, bl, dofVec] = VEM3D_loc_v2(G,f,i);

    S(dofVec,dofVec) = S(dofVec,dofVec) + Sl;
    b(dofVec) = b(dofVec) + bl;

end

[bcDof, bBC] = VEM3D_bc(G,bc);
b(bcDof == 1) = bBC(bcDof == 1);
SBC = spdiags(ones(N,1),0,N,N);
S(bcDof == 1,:) = SBC(bcDof == 1,:);
b(bcDof == 2) = b(bcDof == 2) + bBC(bcDof == 2);