function [Sl]=locS(X)
%--------------------------------------------------------------------------
%   Generates local matrix for diffusion problem
%
%   X:  mx2 matrix containing the m vertices of the polygon.
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
hK = sqrt(2);               %   Element Diameter (NG)

                            %   Monomials. m(i) = m_i
m =      @(X) 1/hK.*[X(:,1)-xK, ...                             %   (1,0)
                     X(:,2)-yK, ...                             %   (0,1)
                     (X(:,1)-xK).^2./hK, ....                   %   (2,0)
                     (X(:,1)-xK).*(X(:,2)-yK)./hK, ....         %   (1,1)
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
                            %   Monomial values at the midpoints.
Xmid = (X(1:end,:) + [X(2:end,:);X(1,:)])/2;
D(n+1:2*n,:) = [ones(n,1), m(Xmid)];
                            %   Integral of monomials over K using
                            %   Gauss-Lobatto quadrature. (NG)
D(NK,:) = [1/36*sum(D(1:n,:)) + 1/9*sum(D(n+1:2*n,:))];
D(NK,1) = D(NK,1) + 4/9;

%%  BUILD MATRIX B.

B=zeros(nk,NK);
                            %   Contribution from boundary integrals.
                            %   Basisfunctions with non-zero degrees of
                            %   freedom on edge j are \phi_j,
                            %   \phi_{mod(j,n)+1} (endpoints) and
                            %   \phi_{j+4} (midpoint). For each edge, 
                            %   \int_\patrial K \nabla m_\alpha \cdot nVec
                            %   is evaluated by Gauss-Lobatto.
for j=1:n                  
                            %   Outward normal of edge j, nVec
    nVec = [X(mod(j,n)+1,2) - X(j,2) , - X(mod(j,n)+1,1) + X(j,1) ]';
                            %   First point of edge j.
    grad_mDOTnVec = grad_m(X(j,:))*nVec;
    B(2:nk,j) = B(2:nk,j) + 1/6.*grad_mDOTnVec;
                            %   Last point of edge j.
    grad_mDOTnVec = grad_m(X(mod(j,n)+1,:))*nVec;
    B(2:nk,mod(j,n)+1) = B(2:nk,mod(j,n)+1) + 1/6.*grad_mDOTnVec;
                            %   Midpoint of edge j.
    grad_mDOTnVec = grad_m(Xmid(j,:))*nVec;
    B(2:nk,j+4) = B(2:nk,j+4) + 2/3.*grad_mDOTnVec;
end
                            %   Contribution from surface integrals.
                            %   \Dela m_\alpha = 2/hK^2 for \alpha = (2,0)
                            %   and (0,2), 0 otherwise. (NG)
B(4,NK) = B(4,NK) - 2/hK^2;
B(6,NK) = B(6,NK) - 2/hK^2;
B(1,NK) = 1;                %   First row (!?). (NG)

%%  CONSTRUCTION OF LOCAL STIFFNESS MATRIX.                              %% 

G = B*D;
PNstar = G\B;               %   \Pi^\Nabla in the monomial basis.
PN = D*PNstar;              %   \Pi^\Nabla in the V^K basis.
Gtilde = zeros(nk,nk);
Gtilde(2:nk,:)=G(2:nk,:);

                            %   Local stiffness matrix.
Sl = PNstar'*Gtilde*PNstar + (eye(NK)-PN)'*(eye(NK)-PN);

D
B
G

end