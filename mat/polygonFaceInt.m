function [IB, ID] = polygonFaceInt(X,hK, Kc, hF, Fc, T,faceNormal,grad_m3D,m2D, int_m3D, PNstar)
%--------------------------------------------------------------------------
%   Evaluates integral of function m : R^2 -> R over a general polygon
%   with vertices X by mapping to reference triangle tRef with vertices
%   {(0,0), (1,0), (0,1)} and using Gauss-Lobatto quadrature. Exact for
%   polynomials of degree <= 2.
%
%   X:  Matrix of vertex coordinates, in counter-clockwise order. Vertex i
%       has coordinates (X(i,1), X(i,2)).
%   m:  Function whose surface integral over the polygn is to evaluated.
%
%--------------------------------------------------------------------------

                            %   Triangulate polygon.
tri = delaunay(X);
nTri = size(tri,1);

                            %   Gauss-Lobatto quadrature point and
                            %   wheights for refenrence triangle.
Xq = [0.0, 0.0; 0.0, 1.0; 0.5, 0.0; 0.5, 0.5; 0.0, 0.5; 0.5, 0.25];
w = [1/36; 1/36; 1/18; 1/18; 1/9; 2/9];

IB = zeros(9,size(PNstar,2));
ID = zeros(1,10);
for t = 1:nTri
                            %   Triangle points
    x1 = X(tri(t,1),:); x2 = X(tri(t,2),:); x3 = X(tri(t,3),:);
                            %   map phi : tRef -> t
    phi = @(x) x*[x1(1)-x3(1)   x1(2)-x3(2);
                  x2(1)-x3(1)   x2(2)-x3(2)] + repmat(x3,size(x,1),1);
                            %   determinant of Dphi, Jacobian of phi.
    detDphi = (x1(1)-x3(1))*(x2(2)-x3(2)) - (x1(2)-x3(2))*(x2(1)-x3(1));
                            %   Map quadrature points to t
    XIq = phi(Xq);
    XK = (XIq*T'+repmat(Fc,6,1)-repmat(Kc,6,1))/hK;
    XF = XIq/hF;
    GRAD_M3D = sum(grad_m3D(XK).*repmat(faceNormal,54,1),2);
    
                            %   Wheights must be scaled by edge lengths.
    GRAD_M3D = bsxfun(@times, reshape(GRAD_M3D,6,9), w);
    GRAD_M3D = GRAD_M3D./hK;
    
    IB = IB + abs(detDphi).*GRAD_M3D'*m2D(XF)*PNstar;  
end

end

%   For each triangle t, evaluate integral.
% I = 0;
% 
%                             %   Triangle points
% x1 = X(tri(:,1),:); x2 = X(tri(:,2),:); x3 = X(tri(:,3),:);
% 
%                             %   map phi : tRef -> t
% phi = zeros(2*nTri,2);
% phi(1:2:2*nTri-1,:) = [x1(:,1)-x3(:,1)   x1(:,2)-x3(:,2)];
% phi(2:2:2*nTri,:)   = [x2(:,1)-x3(:,1)   x2(:,2)-x3(:,2)];
% 
% 
% %     phi = @(x) x*[x1(1)-x3(1)   x1(2)-x3(2);
% %                   x2(1)-x3(1)   x2(2)-x3(2)] + repmat(x3,size(x,1),1);
%                             %   determinant of Dphi, Jacobian of phi.
% detDphi = (x1(:,1)-x3(:,1)).*(x2(:,2)-x3(:,2)) - (x1(:,2)-x3(:,2)).*(x2(:,1)-x3(:,1));
%                         %   Map quadrature points to t
%  
% XIq = Xq*[x1(:,1)-x3(:,1)   x1(:,2)-x3(:,2);
%           x2(:,1)-x3(:,1)   x2(:,2)-x3(:,2)] + repmat(x3,size(Xq,1),1);
% XIqF = XIq*T;
%                         %   Evaluate integral
% % I = I + abs(detDphi)*w*grad_m(XIq)