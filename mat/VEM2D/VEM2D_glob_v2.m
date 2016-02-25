function [A,b] = VEM2D_glob_v2(G,bc,k)

nN = G.nodes.num;
nE = G.faces.num;
nK = G.cells.num;

N = nN + nE*(k-1) + nK*k*(k-1)/2;

if k == 1
    dofPosS = [1, cumsum((diff(G.cells.nodePos')).^2) + 1];
    dofPosb = [1, cumsum(diff(G.cells.nodePos')) + 1];
elseif k == 2
    dofPosS = [1, cumsum((diff(G.cells.nodePos') +  ...
                          diff(G.cells.facePos') + 1).^2) + 1];
    dofPosb = [1, cumsum(diff(G.cells.nodePos') + ...
                         diff(G.cells.facePos') + 1) + 1];
end        

iiS = zeros(1,dofPosS(end)-1);
jjS = zeros(1,dofPosS(end)-1);

sVec = zeros(1,dofPosS(end)-1);
bVec = zeros(1,dofPosb(end)-1);

fprintf('Computing local block matrices ...\n')
tic;

for K = 1:nK
    
    AK = G.cells.AK{K};
    bK = G.cells.bK{K};
    
    dofVec = VEM2D_loc_v2(G,K,k);

    NK = numel(dofVec);
    
    iil = repmat(dofVec', NK, 1);
    jjl = repmat(dofVec , NK, 1);
    jjl = jjl(:);
    
    iiS(dofPosS(K):dofPosS(K+1)-1) = iil;
    jjS(dofPosS(K):dofPosS(K+1)-1) = jjl;
    sVec(dofPosS(K):dofPosS(K+1)-1)= AK(:);

    iib(dofPosb(K):dofPosb(K+1)-1) = dofVec;
    bVec(dofPosb(K):dofPosb(K+1)-1) = bK(:);
    
end

stop = toc;
fprintf('Done in %f seconds.\n\n', stop);

A = sparse(iiS, jjS, sVec, N, N);
b = sparse(iib, ones(1, numel(iib)), bVec, N, 1);

