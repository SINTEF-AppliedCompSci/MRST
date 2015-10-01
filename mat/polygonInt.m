function I = polygonInt(X,m)
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

                            %   For each triangle t, evaluate integral.
I = zeros(1,5);
for t = 1:nTri
                            %   Triangle points
    x1 = X(tri(t,1),:); x2 = X(tri(t,2),:); x3 = X(tri(t,3),:);
                            %   map phi : tRef -> t
    phi = @(x) x*[x1(1)-x3(1)   x1(2)-x3(2);
                  x2(1)-x3(1)   x2(2)-x3(2)] + repmat(x3,size(x,1),1);
                            %   determinant of Dphi, Jacobian of phi.
    detDphi = (x1(1)-x3(1))*(x2(2)-x3(2)) - (x1(2)-x3(2))*(x2(1)-x3(1));
                            %   Gauss-Lobatto quadrature point and
                            %   wheights for refenrence triangle.
    Xq = [0.0, 0.0; 0.0, 1.0; 0.5, 0.0; 0.5, 0.5; 0.0, 0.5; 0.5, 0.25];
    w = [1/36, 1/36, 1/18, 1/18, 1/9, 2/9];
                            %   Map quadrature points to t
    XIq = phi(Xq);
                            %   Evaluate integral
    I = I + abs(detDphi)*w*m(XIq);
end

end