function [Sl, bl] = locS(X, Xmid, edgeLengths, edgeNormals, onBoundary, vol, f,hK)
%--------------------------------------------------------------------------
%   Generates local matrix for diffusion problem
%
%   X:  nx2 matrix containing the n vertices of the polygon.
%       X(i,:) = (x,y) coordinates of vertex i. Vertices ordered counter
%       clockwise.
%   Sl: NK x NK local stifness matrix.
%--------------------------------------------------------------------------

%%  DIMENSIONS, BARICENTER, MONOMIALS.                                   %% 

n = size(X,1);              %   Number of vertices.
k = 2;                      %   Method order.
nk = 0.5*(k+1)*(k+2);       %   Dimension of polynomial space.
NK = n*k + 0.5*k*(k-1);     %   Dimension of V^K (degrees of freedom).
[~,BB]=baric(X);            %   Baricenter of K
xK = BB(1); yK = BB(2);        
%hK = 5;               %   Element Diameter (NG)

                            %   Monomials. m(i) = m_i
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
%D(NK,:) = [1/36*sum(D(1:n,:)) + 1/9*sum(D(n+1:2*n,:))];
%D(NK,1) = D(NK,1) + 4/9;

%%  BUILD MATRIX B.

B=zeros(nk,NK);
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
                            %   First point of edge j.
        grad_mDOTnVec = grad_m(X(j,:))*nVec;
        B(2:nk,j) = B(2:nk,j) + edgeLengths(j)*1/6.*grad_mDOTnVec;
                            %   Last point of edge j.
        grad_mDOTnVec = grad_m(X(mod(j,n)+1,:))*nVec;
        B(2:nk,mod(j,n)+1) = B(2:nk,mod(j,n)+1) + edgeLengths(j)*1/6.*grad_mDOTnVec;
                            %   Midpoint of edge j.
        grad_mDOTnVec = grad_m(Xmid(j,:))*nVec;
        B(2:nk,j+n) = B(2:nk,j+n) + edgeLengths(j)*2/3.*grad_mDOTnVec;

end
                            %   Contribution from surface integrals.
                            %   \Delta m_\alpha = 2/hK^2 for \alpha = (2,0)
                            %   and (0,2), 0 otherwise.
B(4,NK) = B(4,NK) - vol*2/hK^2;
B(6,NK) = B(6,NK) - vol*2/hK^2;
B(1,NK) = 1;                %   First row (!?).

%%  CONSTRUCTION OF LOCAL STIFFNESS MATRIX.                              %% 

G = B*D;
PNstar = G\B;               %   \Pi^\Nabla in the monomial basis.
PN = D*PNstar;              %   \Pi^\Nabla in the V^K basis.
Gtilde = [zeros(1,nk) ; G(2:nk,:)];

                            %   Local stiffness matrix.
Sl = PNstar'*Gtilde*PNstar + (eye(NK)-PN)'*(eye(NK)-PN);

bl = zeros(NK,1);
% for j = 1:n
%     if ~onBoundary(j)
%         bl(j) = bl(j) + 1/6*f(X(j,:));
%         bl(mod(j,n)+1) = bl(mod(j,n)+1) + 1/6*f(X(mod(j,n)+1,:));
%         bl(j+4) = bl(j+4) + 2/3.*f(Xmid(j,:));
%     end
% end
        
end