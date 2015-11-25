function [Sl, bl, dofVec, hK] = VEM2D_loc(G, K, f)
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
[nodes, X] = nodeData(G,K);
                            %   Cell edges and edge midpoint coordinates.
[edges, Xmid, edgeNormals] = faceData(G,K);
                            %   Baricenter of K.
Xc = G.cells.centroids(K,:); xK = Xc(1); yK = Xc(2);
hK = cellDiameter(X);       %   Element diameter.
vol = G.cells.volumes(K);   %   Cell volume.
n = size(X,1);              %   Number of vertices.
k = 2;                      %   Method order.
nk = 0.5*(k+1)*(k+2);       %   Dimension of polynomial space.
NK = n*k + 0.5*k*(k-1);     %   Local nomber of dofs.

                            %   Monomials. m(i) = m_{i+1}.
m =      @(X) [(X(:,1)-xK)./hK                , ...          %   (1,0)
               (X(:,2)-yK)./hK                , ...          %   (0,1)
               (X(:,1)-xK).^2./hK^2           , ...          %   (2,0)
               (X(:,1)-xK).*(X(:,2)-yK)./hK^2 , ...          %   (1,1)
               (X(:,2)-yK).^2./hK^2            ];            %   (0,2)
                 
                            %   Gradients of monomials.
                            %   grad_m(i,:) = \nabla m_{i+1}.
grad_m = @(X) [1/hK                , 0                  ;    %   (1,0)
               0                   , 1/hK               ;    %   (0,1)
               (X(:,1)-xK).*2/hK^2 , 0                  ;    %   (2,0)
               (X(:,2)-yK)./hK^2   , (X(:,1)-xK)./hK^2  ;    %   (1,1)
               0                   , (X(:,2)-yK).*2/hK^2 ];  %   (0,2)
           
m_int = @(X)  [(X(:,1)-xK).^2./(2*hK) , ...
               (X(:,1)-xK).*(X(:,2)-yK)./hK , ...
               (X(:,1)-xK).^3./(3*hK^2) , ...
               (X(:,1)-xK).^2.*(X(:,2)-yK)./(2*hK^2) , ...
               (X(:,1)-xK).*(X(:,2)-yK).^2./(hK^2)];
                
%%  BUILD MATRIX D                                                       %%

                            %   Integral of monomials over K using
                            %   Gauss-Lobatto quadrature.
% I = [vol, polygonInt(X,m)];
I = [vol, evaluateMonomialIntegralV2(edgeNormals, X, Xmid, m_int)];
                            %   D(i,j) = \chi^i(m_j(X_i))
D = [ones(NK-1,1), m([X;Xmid]); I./vol];

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
                            %   Area wheighted outward normal of edge j.
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

                            %   Matrix of integrals over K of all
                            %   combinations of linear monomials.
H = [I(1:3)        ;
     I(2)   I(4:5) ;
     I(3)   I(5:6)];
% H = H + tril(H',-1);s
                            %   \Pi^\Nabla in the monomial basis
                            %   \mathcal{M}_1.
PNstar = M(1:3,1:3)\B(1:3,:);
                            %   Local load term.
bl = PNstar'*H*PNstar*f([X; Xmid; Xc]);

%%  CONSTRUCT LOCAL TO GLOBAL MAP. S(dofVec,dofVec) = Sl.                %%

                            %   Dofs are ordered in the following way:
                            %   values at nodes, values at edge midpoints,
                            %   values of moments on faces.

Nn = G.nodes.num;           %   Total number of nodes.
Ne = G.faces.num;           %   Total number of edges.

nodes = fixDim(nodes);      %   Fix dimensions.
edges = fixDim(edges);
                            %   Local to global map,
                            %   S(dofVec, dofVec) = Sl.
dofVec = [nodes, Nn + edges, Nn + Ne + K];

end