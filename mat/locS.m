function [lS]=locS(X)
%--------------------------------------------------------------------------
%   Generates local matrix for diffusion problem
%
%   X:  mx2 matrix containing the m vertices of the polygon.
%       X(i,:) = (x,y) coordinates of vertex i. Vertices ordered counter
%       clockwise.
%   Sl: m x m local stifness matrix.
%--------------------------------------------------------------------------

                            %   Number of vertices.                       %
n = size(X,1);
k = 2;
nk = 0.5*(k+1)*(k+2);
NK = n*k + 0.5*k*(k-1);
alf=1;
[~,BB]=baric(X);            %   Baricenter coordinates.                  %
xV = X(:,1); yV = X(:,2);
xK = BB(1); yK = BB(2);

% -------
% BASIS for the VEM space is value at first
% vertex, then second vertex, etc ... 
% BASIS for the polynomial space is 
% {1,x,y} with coordinates x,y scaled with the baricenter
                            
                            %   Build matrix D.                           %
hK = sqrt(2);               %   Element Diameter
                            %   Monomials and their gradients and         %
                            %   Laplacians.                               %
m1 = @(x,y) 1*(x+1)./(x+1);                      %   m_(0,0)
m2 = @(x,y) (x-xK)./hK;             %   m_(1,0)
m3 = @(x,y) (y-yK)./hK;             %   m_(0,1)
m4 = @(x,y) (x-xK).^2/hK^2;         %   m_(2,0)
m5 = @(x,y) (x-xK).*(y-yK)./hK^2;   %   m_(1,1)
m6 = @(x,y) (y-yK).^2./hK^2;        %   m_(0,2)

grad_m2 = @(x,y) 1/hK*[1,0];        %   \nabla m_(1,0)
grad_m3 = @(x,y) 1/hK*[0,1];        %   \nabla m_(0,1)
grad_m4 = @(x,y) 1/hK^2*[x-xK,0];   %   \nabla m_(2,0)
grad_m5 = @(x,y) 1/hK^2*[1,1];      %   \nabla m_(1,1)
grad_m6 = @(x,y) 1/hK^2*[0,y-yK];   %   \nabla m_(0,2)

malpha = (X - ones(n,1)*BB)./hK;
D = zeros(NK, nk);
m2(xV,yV)
                            %   Monomial values at the vertices.          %
D(1:n,:) = [m1(xV,yV), m2(xV,yV), m3(xV,yV), ...
            m4(xV,yV), m5(xV,yV), m6(xV,yV)];
                            %   Monomial values at the midpoints.         %
mid = (X(1:end,:) + [X(2:end,:);X(1,:)])/2;
xMid = mid(:,1); yMid = mid(:,2);
D(n+1:2*n,:) = [m1(xMid,yMid), m2(xMid,yMid), m3(xMid,yMid), ...
                m4(xMid,yMid), m5(xMid,yMid), m6(xMid,yMid)];
                            %   Integral of monomials over the element    %
                            %   using Gauss-Lobatto quadrature.           %
D(NK,:) = [1/36*sum(D(1:n,:)) + 1/9*sum(D(n+1:2*n,:))];
D(NK,1) = D(NK,1) + 4/9


                            %   Build matrix B (initially as transpose).  %
B=zeros(NK,nk);
for j=1:nk                  %   Add contribution of each edge             %
                            %   Kernel part (first column).
                            %   See projection below.                     %
   B(j,1) = (1/n); 
                            %   Non kernel part (second and third columns)%
                            %   Outward normal.                           %
   n = [ X(mod(j,m)+1,2) - X(j,2) , - X(mod(j,m)+1,1) + X(j,1) ];
   B(j,2:3) = B(j,2:3) + (1/2)*n;
   B(mod(j,m)+1,2:3) = B(mod(j,m)+1,2:3) + (1/2)*n; 
end
B=B'; % transpose 

% --- build projectors
G = B*D;          % G matrix from Hitchhikers paper
PNs = inv(G)*B;   % PiNabla star projector (polynomial basis)
PN = D*PNs;      % PiNabla projector (Vh basis)  
%  --- build local matrix ---
Gt(1,:)=zeros(1,3);   % matrix G tilde (null on the kernel)
Gt(2:3,:)=G(2:3,:);

M = PNs'*Gt*PNs;  % consistent part
Sl = M + alf*(trace(M)/2)*(eye(m)-PN)'*(eye(m)-PN);


% -----------------------
% KERNEL PROJECTION
% -----------------------
% Let the vertex based scalar product 
% <v,w> = (1/m) sum_{i=1}^m v(vertex_i)\cdot w(vertex_i) .
% Then, the projector on the kernel is
% P(v) := <1,v>
% -----------------------