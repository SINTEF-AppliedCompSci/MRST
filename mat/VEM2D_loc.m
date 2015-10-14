function [Sl, bl, dofVec] = VEM2D_loc(G, K, f)
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
%   Output:
%
%   Sl:     Local stiffness matrix for cell K. dim(Sl) = NK x NK, where
%           NK = n*k + 0.5*k*(k-1), n is the number of vertices of K,
%           and k = 2 is the order of the method.
%   
%   dofVec: Map from local to global dofs. S(dofVec, dofVec) = Sl, where
%           S is the global stiffness matrix.
%   
%--------------------------------------------------------------------------

%%  CELL DATA                                                            %%

                            %   Cell nodes and node coordinates.
nodeNum = G.cells.nodePos(K):G.cells.nodePos(K+1)-1;
cellNodes = G.cells.nodes(nodeNum);
X = G.nodes.coords(cellNodes,:);
                            %   Cell edges and edge midpoint coordinates.
faceNum = G.cells.facePos(K):G.cells.facePos(K+1)-1;
cellFaces = G.cells.faces(faceNum);
Xmid = G.faces.centroids(cellFaces,:);
                            %   Area wheighted edge normals.
edgeNormals = G.faces.normals(cellFaces,:);
neighbors = G.faces.neighbors(cellFaces,:);
                            %   Fix orientation of edge normals.
m = (-ones(length(edgeNormals),1)).^(neighbors(:,1) ~= K);
edgeNormals = [m,m].*edgeNormals;

vol = G.cells.volumes(K);   %   Cell volume.
n = size(X,1);              %   Number of vertices.
k = 2;                      %   Method order.
nk = 0.5*(k+1)*(k+2);       %   Dimension of polynomial space.
NK = n*k + 0.5*k*(k-1);     %   Dimension of V^K (degrees of freedom).
[~,XB]=baric(X);            %   Baricenter of K.
xK = XB(1); yK = XB(2);

                            %   Element diameter.
hK = 0;
for i = 1:n
    hK = max(norm(repmat(X(i,:),n,1)-X),hK);
end

                            %   Monomials. m(i) = m_i.
m =      @(X) 1/hK.*[X(:,1)-xK, ...                             %   (1,0)
                     X(:,2)-yK, ...                             %   (0,1)
                     (X(:,1)-xK).^2./hK, ....                   %   (2,0)
                     (X(:,1)-xK).*(X(:,2)-yK)./hK, ...          %   (1,1)
                     (X(:,2)-yK).^2./hK];                       %   (0,2)
                 
                            %   Gradients of monomials.
                            %   grad_m(i,:) = \nabla m_i.
grad_m = @(X) 1/hK.*[1                 , 0;                     %   (1,0)
                     0                 , 1;                     %   (0,1)
                     2/hK.*(X(:,1)-xK) , 0;                     %   (2,0)
                     1/hK.*(X(:,2)-yK) , 1/hK.*(X(:,1)-xK);     %   (1,1)
                     0                 , 2/hK.*(X(:,2)-yK)   ]; %   (0,2)

%%  BUILD MATRIX D                                                       %%

D = zeros(NK, nk);
                            %   Monomial values at the vertices.
D(1:n,:) = [ones(n,1), m(X)];
                            %   Monomial values at the edge midpoints.
D(n+1:2*n,:) = [ones(n,1), m(Xmid)];
                            %   Integral of monomials over K using
                            %   Gauss-Lobatto quadrature.
D(NK,:) = [1, polygonInt(X,m)./vol];

%%  BUILD MATRIX B.

B=zeros(nk,NK);
B(1,NK) = 1;                %   First row.

                            %   Contribution from boundary integrals.
                            %   Basisfunctions with non-zero degrees of
                            %   freedom on edge j are \phi_j,
                            %   \phi_{mod(j,n)+1} (endpoints) and
                            %   \phi_{j+n} (midpoint). For each edge, 
                            %   \int_\patrial K \nabla m_\alpha \cdot nVec
                            %   is evaluated by Gauss-Lobatto quadrature.
for j=1:n
                            %   Outward normal of edge j, nVec
    nVec = edgeNormals(j,:)';
    B(2:nk,[j, mod(j,n)+1, j + n]) = B(2:nk,[j, mod(j,n)+1, j + n]) + ...
           [1/6.*grad_m(X(j,:))*nVec          , ...    %   First point.
            1/6.*grad_m(X(mod(j,n)+1,:))*nVec , ...    %   Last point.
            2/3*grad_m(Xmid(j,:))*nVec             ];  %   Midpoint.
end
                            %   Contribution from surface integrals.
                            %   \Delta m_\alpha = 2/hK^2 for \alpha = (2,0)
                            %   and (0,2), 0 otherwise.
B([4, 6], NK) = B([4, 6], NK) - vol*2/hK^2;

%%  CONSTRUCTION OF LOCAL STIFFNESS MATRIX.                              %% 

M = B*D;
PNstar = M\B;               %   \Pi^\Nabla in the monomial basis.
PN = D*PNstar;              %   \Pi^\Nabla in the V^K basis.
Mtilde = [zeros(1,nk) ; M(2:nk,:)];

                            %   Local stiffness matrix.
Sl = PNstar'*Mtilde*PNstar + (eye(NK)-PN)'*(eye(NK)-PN);

%%  LOCAL LOAD TERM.                                                     %%

bl = zeros(NK, 1);
k = 1;                       %   Method order.
nk1 = 0.5*(k+1)*(k+2);       %   Dimension of polynomial space.
NK1 = NK; %n*k + 0.5*k*(k-1);     %   Dimension of V^K (degrees of freedom).
% B = B(1:nk1, 1:NK1); B(1,:) = 1/NK1.*ones(1,NK1);
% D = D(1:NK1, 1:nk1);
% M = B*D;
% PNstar = M\B;
bl(NK) = vol*PNstar(nk,:)*f([X; Xmid; XB] );


%%  CONSTRUCT LOCAL TO GLOBAL MAP. S(dofVec,dofVec) = Sl.                %%

                            %   Dofs are ordered in the following way:
                            %   values at nodes, values at edge midpoints,
                            %   values of moments on faces.

Nn = G.nodes.num;           %   Total number of nodes.
Ne = G.faces.num;           %   Total number of edges.

                            %   Fix dimensions.
if size(cellNodes, 2) == 1
    cellNodes = cellNodes';
end
if size(cellFaces, 2) == 1
    cellFaces = cellFaces';
end
dofVec = [cellNodes, Nn + cellFaces, Nn + Ne + K];


end