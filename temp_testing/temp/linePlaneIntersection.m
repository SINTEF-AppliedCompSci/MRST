function [Int,type] = linePlaneIntersection(n,p,pline)
Int = NaN(1,3);

p1 = pline(1,:);
p2 = pline(2,:);

l_vec = p1-p2;
pl_vec = p1-p;

d = dot(n,l_vec);
N = -dot(n,pl_vec);

if abs(d) < eps*100 % The segment is parallel to plane
    if N == 0 % The segment lies in plane
        type = 2;
        Int = pline;
        return
    else % no intersection
        type = 0; 
        return
    end
else
    % compute the intersection parameter
    sI = N/d;
    Int = round(p1 + sI.*l_vec, 14);
    if (sI < 0 || sI > 1) % Int outside line-segment limits.
        type = 3; 
    else
        type = 1;
    end
end
return