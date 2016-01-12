function [Sl, bl, dofVec, hK] = VEM3D_loc(G, K, f)
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
%           G = globalEdgeData(G);
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

                            %   Cell nodes and node coordinates.
[vertices, X] = nodeData(G,K);
                            %   Cell faces, face midpoints and normals.
[faces, Fc, faceNormals] = faceData(G,K);
                            %   Cell edges, edge midpoints and normals.
[edges, Ec, edgeNormals] = edgeData(G,faces);

% nodeNum = mcolon(G.edges.nodePos(edges),G.edges.nodePos(edges+1)-ones(size(edges,1),1))
% nodes = G.edges.nodes(nodeNum)
% plot3(G.nodes.coords(nodes,1), G.nodes.coords(nodes,2), G.nodes.coords(nodes,3), '*')

                            %   Baricenter of K.
Kc = G.cells.centroids(K,:); xK = Kc(1); yK = Kc(2); zK = Kc(3);
hK = cellDiameter(X);       %   Element diameter.
vol = G.cells.volumes(K);   %   Element volume.

nV = size(X,1);             %   Number of vertices.
nE = size(edges,1);         %   Number of edges.
nF = size(faces,1);         %   Number of faces.

k = 2;                      %   Method order.
nk = (k+1)*(k+2)*(k+3)/6;   %   Dimension of polynomial space.
                            %   Local nomber of dofs.
NK = nV + nE*(k-1) + nF*k*(k-1)/2 + k*(k^2-1);

%%  BUILD MATRIX Df
    
           
           
%%  BUILD MATRIX D                                                       %%



                            %   Integral of monomials over K using
                            %   Gauss-Lobatto quadrature.
I = [vol, polygonInt2(edgeNormals, X, Xmid, m_int)];
                            %   D(i,j) = \chi^i(m_j(X_i))
D = [ones(NK-1,1), m([X;Xmid]); I./vol];
%%  BUILD MATRIX B.

B=zeros(nk,NK);
B(1,NK) = 1;                %   First row.

                            %   Contribution from boundary integrals.
                            %   Basisfunctions with non-zero degrees of
                            %   freedom on edge j are \phi_j,
                            %   \phi_{mod(j,n)+1} (endpoints) and (nodes);
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

                            %   Local to global map,
                            %   S(dofVec, dofVec) = Sl.
dofVec = [nodes, Nn + edges, Nn + Ne + K];

end