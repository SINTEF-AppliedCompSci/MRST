function I = evaluateMonomialIntegralV2(X,m)

n = size(X,1);
                            %   Triangulate polygon
tri = delaunay(X);
nTri = size(tri,1);

I = 0;
for t = 1:nTri
    pts = X(tri(t,:),:)
    [~,i] = sort(pts,1);
    pts = pts(i(:,1),:);
    tmpy = @(x) (x-pts(1,1)).*(pts(3,2)-pts(1,2))/(pts(3,1)-pts(1,1)) + pts(1,2);
    a1 = (pts(2,1)-pts(1,1));
    tmpy1 = @(x) (x-pts(1,1)).*(pts(2,2)-pts(1,2))/(pts(2,1)-pts(1,1)) + pts(1,2);
    a2 = (pts(3,1)-pts(2,1));
    tmpy2 = @(x) (x-pts(2,1)).*(pts(3,2)-pts(2,2))/(pts(3,1)-pts(2,1)) + pts(2,2);
    if pts(2,2) > tmpy(pts(2,1))
        if a1 > eps
            I = I + quad2d(m,pts(1,1),pts(2,1),tmpy, tmpy1);
        end
        if a2 > eps
            I = I + quad2d(m,pts(2,1),pts(3,1),tmpy, tmpy2);
        end 
    else
        if a1 > eps
            I = I + quad2d(m,pts(1,1),pts(2,1),tmpy1, tmpy);
        end
        if a2 > eps
            I = I + quad2d(m,pts(2,1),pts(3,1),tmpy2, tmpy);
        end
    end
end