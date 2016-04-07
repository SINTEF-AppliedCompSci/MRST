function [AK, bK, dofVec] = VEM2D_loc_v3(G, f, K, k, alpha, rate)
%--------------------------------------------------------------------------
%   Generates local stffness matrix for the virtual element method  for
%   cell K of grid G for diffusion problem:
%
%           -\delta u = f, x \in \Omega
%                   u = g, x \in \partial \Omega
%
%   Input:
%
%   G:      2D MRST grid. Cells can be any kind of polygn. the function
%           assumes the following functions has been called for the grid:
%
%           G = mrstGridWithFullMappings(G);
%           G = computeGeometry(G);
%
%   K:      Cell number in grid G, i.e. G.cells(K).
%
%   f:      Source term.
%
%   Output:
%
%   Sl:     Local stiffness matrix for cell K. dim(Sl) = NK x NK, where
%           NK = n*k + 0.5*k*(k-1), n is the number of vertices of K,
%           and k = 2 is the order of the method.
%
%   bl:     Local load vector for cell K. dim(bl) = 1 x NK.
%   
%   dofVec: Map from local to global dofs. S(dofVec, dofVec) = Sl, where
%           S is the global stiffness matrix.
% 
%--------------------------------------------------------------------------

%%  CELL DATA                                                            %%

[m, grad_m, int_m] = retrieveMonomials(k);

Kc = G.cells.centroids(K,:);
hK = G.cells.diameters(K);
aK = G.cells.volumes(K);

edgeNum = G.cells.facePos(K):G.cells.facePos(K+1)-1;
edges   = G.cells.faces(edgeNum);
if size(edges,1) == 1;
    edges = edges';
end
nE          = size(edges,1);
Ec          = G.faces.centroids(edges,:);
hE          = G.faces.areas(edges,:);
edgeNormals = G.faces.normals(edges,:);
edgeSign    = (-ones(nE,1)).^(G.faces.neighbors(edges,1) ~= K); 
edgeNormals = bsxfun(@times, edgeNormals, edgeSign);

nodeNum = mcolon(G.faces.nodePos(edges),G.faces.nodePos(edges+1)-1);
nodes   = G.faces.nodes(nodeNum);
if size(nodes,1) == 1
    nodes = nodes';
end
nodes   = reshape(nodes,2,[])';
nN      = size(nodes,1);
nodes(edgeSign == -1,:) = nodes(edgeSign == -1,2:-1:1);
nodes   = nodes(:,1);

nk = (k+1)*(k+2)/2;
nkk = k*(k+1)/2;
NK = nN + nE*(k-1) + k*(k-1)/2;

X = [G.nodes.coords(nodes,:); Ec];

Xmon = bsxfun(@minus,X,Kc)/hK;                     

%%  CALCULATE B AND D MATRICES                                           %%

if k == 1
    
    D = m(Xmon(1:nN,:));
    
    XintB = repmat(Xmon(nN+1:end,:),2,1);
    edgeNormals = [edgeNormals;...
                   [edgeNormals(nE,:);edgeNormals(1:nE-1,:)]];
    
    
    intB = .5*sum(grad_m(XintB).*repmat(edgeNormals,2,1),2);
    intB = reshape(intB,2*nE,2);
    intB = (intB(1:nE,:) + intB(nE+1:2*nE,:))/hK;
    B = [ones(1,NK)/NK; intB'];
    
    H = aK;
    
    if isa(f,'function_handle')
        fHat = f(X(1:nN,:)) + rate/aK;
    else
        fHat = f*ones(size(nN,1)) + rate/aK;
    end
    
    dofVec = nodes
    
elseif k == 2
    
    intD = bsxfun(@times                               , ...
                  (int_m(Xmon(1:nN,:))                          ...
                   + [int_m(Xmon(2:nN,:));int_m(Xmon(1,:))])/6  ...
                   +  int_m(Xmon(nN+1:end,:))*2/3              , ...
                   edgeNormals(:,1));
    intD = sum(intD,1)*hK/aK;
    D = [m(Xmon); intD];

    XintB = [Xmon(1:nN,:); Xmon(1:nN,:); Xmon(nN+1:end,:)];
    edgeNormals = [edgeNormals;...
                   [edgeNormals(nE,:);edgeNormals(1:nE-1,:)]; ...
                   edgeNormals];
    
    intB = sum(grad_m(XintB).*repmat(edgeNormals,5,1),2);
    intB = reshape(intB,3*nN,5);
    intB = [(intB(1:nN,:) + intB(nN+1:2*nN,:))/6; ...
             intB(2*nN+1:end,:)*2/3]/hK;
         
    B = zeros(nk, NK);
    B(1,NK) = 1;
    B(2:nk, 1:NK-1) = intB';
    B([4, 6], NK) = B([4, 6], NK) - 2*aK/hK^2;

    H = zeros(nkk, nkk);
    H(1,:) = intD(:, [1,2,3])*aK;
    H(2,:) = intD(:, [2,4,5])*aK;
    H(3,:) = intD(:, [3,5,6])*aK;
    
    if isa(f,'function_handle')
        fInt = polygonInt_v2(G, K, f, k+1)/aK;
        fHat = [f(X); fInt] + rate/aK;
    else
        fHat = f*ones(2*nN+1,1) + rate/aK;
    end
    
    dofVec = [nodes', edges' + G.nodes.num, K + G.nodes.num + G.faces.num];
    
end

M = B*D;
PNstar = M\B;
PN = D*PNstar;
% Q = sqrt(9/(4*.5*hK^2))*orth(eye(NK) - PN);
% P = Q'*Q;
Mtilde = [zeros(1,nk) ; M(2:nk,:)];
% AK = PNstar'*Mtilde*PNstar ...
%            + alpha(K)*(eye(NK)-PN)'*(Q/P)*(P\Q')*(eye(NK)-PN);
AK = PNstar'*Mtilde*PNstar ...
           + alpha(K)*(eye(NK)-PN)'*(eye(NK)-PN);
PNstar = M(1:nkk,1:nkk)\B(1:nkk,:);
bK = PNstar'*H*PNstar*fHat;

end