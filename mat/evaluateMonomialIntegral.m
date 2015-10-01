function I = evaluateMonomialIntegral(X,m)

n = size(X,1);
                            %   Triangulate polygon
tri = delaunay(X);
nTri = size(tri,1);

I = 0;
for t = 1:nTri
    x1 = X(tri(t,1),:);
    x2 = X(tri(t,2),:);
    x3 = X(tri(t,3),:);
    detJ = 1/((x1(1)-x3(1))*(x2(2)-x3(2)) - (x1(2)-x3(2))*(x2(1)-x3(1)));
    A = detJ*[(x2(2)-x3(2))  -(x2(1)-x3(1)); -(x1(2)-x3(2)) x1(1)-x3(1)];
    I = I + detJ*(1/36*(m([0,0]) + m([0,1])) + 1/18*(m([0.5,0]) + ...
                m([0.5,0.5])) + 1/9*m([0,0.5]) + 2/9*m([0.5, 0.25]));
end

m = @(x,y) (x- 1/2).*(y-1/2);
1/36*(m(0,0) + m(0,1)) + 1/18*(m(0.5,0) + ...
                m(0.5,0.5)) + 1/9*m(0,0.5) + 2/9*m(0.5, 0.25)
quad2d(m, 0,1,0,1)
end