function [S,b] = VEM3D_glob(G,f,bc)

nN = G.nodes.num;
nE = G.edges.num;
nF = G.faces.num;
nK = G.cells.num;
k = 2;

N = nN + nE*(k-1) + nF*k*(k-1)/2 + nK*k*(k^2-1)/6;

dofPosS = [1, cumsum((diff(G.cells.nodePos') + ...
                      diff(G.cells.edgePos') + ...
                      diff(G.cells.facePos') + 1).^2) + 1];
dofPosb = [1, cumsum(diff(G.cells.nodePos') + ...
                     diff(G.cells.edgePos') + ...
                     diff(G.cells.facePos') + 1) + 1];
                 
step = floor(nK/10);

% itS = 1; itb = 1;

iiS = zeros(1,dofPosS(end)-1);
jjS = zeros(1,dofPosS(end)-1);

sVec = zeros(1,dofPosS(end)-1);
bVec = zeros(1,dofPosb(end)-1);

for i = 1:nK

    if rem(i,step) == 0
        fprintf('Calculating local block matrix for cell %d ...\n', i);
    end
    
    [Sl, bl, dofVec] = VEM3D_loc(G,f,i);

    NK = numel(dofVec);
    
    iil = repmat(dofVec', NK, 1);
    jjl = repmat(dofVec , NK, 1);
    jjl = jjl(:);
    
    iiS(dofPosS(i):dofPosS(i+1)-1) = iil;
    jjS(dofPosS(i):dofPosS(i+1)-1) = jjl;
    sVec(dofPosS(i):dofPosS(i+1)-1)= Sl(:);

    iib(dofPosb(i):dofPosb(i+1)-1) = dofVec;
    bVec(dofPosb(i):dofPosb(i+1)-1) = bl(:);
    
end

S = sparse(iiS, jjS, sVec, N, N);
b = sparse(iib, ones(1, numel(iib)), bVec);

[bcDof, bBC] = VEM3D_bc(G,bc);
b(bcDof == 1) = bBC(bcDof == 1);
SBC = spdiags(ones(N,1),0,N,N);
S(bcDof == 1,:) = SBC(bcDof == 1,:);
b(bcDof == 2) = b(bcDof == 2) + bBC(bcDof == 2);