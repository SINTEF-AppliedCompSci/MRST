function cent=findtrisubsegmentcentroid(points1,points2)
% points1 (one or two points) and points2 (two points) are the leftover
% triangle points and the cutting line of a triangular-cell sub-segment.
%
% REVISED: same non-orthonormal-frame defect as convexpolygoncentroid; the
% in-plane frame is made orthonormal using the polygon normal so the centroid
% is order-independent.

points=[points1;points2];
numpoints=size(points,1);

if numpoints<=3
    cent=findtrianglecentroid(points);
    return;
else
    % ORTHONORMAL in-plane frame
    midpoint=sum(points,1)/numpoints;
    normal=cross(points(1,:)-points(3,:),points(2,:)-points(3,:));
    normal=normal/norm(normal);
    axis_x=(points(1,:)-points(3,:)); axis_x=axis_x/norm(axis_x);
    axis_y=cross(normal,axis_x);      axis_y=axis_y/norm(axis_y);

    inplanepointsx=zeros(numpoints,1);inplanepointsy=zeros(numpoints,1);
    for i=1:numpoints
        inplanepointsx(i)=dot(points(i,:)-midpoint,axis_x);
        inplanepointsy(i)=dot(points(i,:)-midpoint,axis_y);
    end

    % choose the sequential arrangement (max area), as in the original
    x1=inplanepointsx([3,4,1,2]); y1=inplanepointsy([3,4,1,2]); p=[2,3,4,1];
    area1=abs((1/2)*sum((x1(p)+x1).*(y1(p)-y1)));
    x2=inplanepointsx([3,4,2,1]); y2=inplanepointsy([3,4,2,1]);
    area2=abs((1/2)*sum((x2(p)+x2).*(y2(p)-y2)));
    if area1>area2
        inplanepointsx=x1; inplanepointsy=y1;
    else
        inplanepointsx=x2; inplanepointsy=y2;
    end

    x=[inplanepointsx;inplanepointsx(1)];
    y=[inplanepointsy;inplanepointsy(1)];

    cx=0;cy=0;area=0;
    for i=1:numpoints
        mult=x(i)*y(i+1)-x(i+1)*y(i);
        cx=cx+(x(i)+x(i+1))*mult;
        cy=cy+(y(i)+y(i+1))*mult;
        area=area+0.5*mult;
    end
    cx=cx/(6*area);
    cy=cy/(6*area);

    cent=cx*axis_x+cy*axis_y+midpoint;
end

end

function tricent=findtrianglecentroid(points)
tricent=sum(points,1)/size(points,1);
end
