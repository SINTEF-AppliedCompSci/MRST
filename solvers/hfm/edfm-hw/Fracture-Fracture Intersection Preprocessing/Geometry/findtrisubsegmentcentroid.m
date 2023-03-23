function cent=findtrisubsegmentcentroid(points1,points2)
% points is an Nx3 matrix with each row representing a vertex and
% containing the xyz coordinates of that vertex. We consider two cases:
% triangle and quadrilaterals.
%
% points1 will contain either one or two points, listed in each row with
% xyz coordinates. points2 will contain two points, listed in each row with
% xyz coordinates. points1 are connected by an edge and so are points2. We
% exploit this property to figure out the sequential arrangement of
% the polygon vertices along the perimeter.
%
% This is mainly for the purpose of triangular cell subsegments. Points2 is
% the line cutting through the triangle. Points1 are the leftover triangle
% points.


points=[points1;points2];

numpoints=size(points,1);

if numpoints<=3
    cent=findtrianglecentroid(points);
    return;
else
    % determine in-plane x-y coordinates of points
    midpoint=sum(points,1)/numpoints;
    axis_x=(points(1,:)-points(3,:))/norm(points(1,:)-points(3,:));
    axis_y=(points(2,:)-points(3,:))/norm(points(2,:)-points(3,:)); %intermediate variable
    axis_z=cross(axis_x,axis_y);
    axis_y=cross(axis_z,axis_x);

    % generate inplanepoints which contains the in-plane x-y coordinates for
    % the 3D points in 'points'. The sequence in inplanepoints corresponds to
    % that in points at the moment.
    inplanepointsx=zeros(numpoints,1);inplanepointsy=zeros(numpoints,1);
    for i=1:numpoints
        inplanepointsx(i)=dot(points(i,:)-midpoint,axis_x);
        inplanepointsy(i)=dot(points(i,:)-midpoint,axis_y);
    end
    
    % test two arrangements here to see which one makes sense. We exploit
    % the fact that we already know two edges of the quadrilateral. So we
    % only need to establish the two other edges. If the two other edges
    % intersect, then the area calculated using the following algorithm
    % will be the difference between the two resulting triangles. If the
    % edges do not intersect, we will obtain a larger area.
    
    % arrangement 1
    % points2(1,:), points2(2,:), points1(1,:), points1(2,:)
    x1=inplanepointsx([3,4,1,2]);
    y1=inplanepointsy([3,4,1,2]);
    p=[2,3,4,1];
    area1=abs((1/2)*sum((x1(p)+x1).*(y1(p)-y1)));
    
    % arrangement 2
    % points2(1,:), points2(2,:), points1(2,:), points1(1,:)
    x2=inplanepointsx([3,4,2,1]);
    y2=inplanepointsy([3,4,2,1]);
    area2=abs((1/2)*sum((x2(p)+x2).*(y2(p)-y2)));
    
    % max area corresponds to a sequential arrangement
    if area1>area2
        inplanepointsx=x1;
        inplanepointsy=y1;
    else
        inplanepointsx=x2;
        inplanepointsy=y2;
    end
    
    x=[inplanepointsx;inplanepointsx(1)];
    y=[inplanepointsy;inplanepointsy(1)];
    
    %generate in-plane centroid location
    cx=0;cy=0;area=0;
    for i=1:numpoints
        mult=x(i)*y(i+1)-x(i+1)*y(i);
        cx=cx+(x(i)+x(i+1))*mult;
        cy=cy+(y(i)+y(i+1))*mult;
        area=area+0.5*mult;
    end
    cx=cx/(6*area);
    cy=cy/(6*area);
    
    % convert in plane coordinates to original coordinates
    cent=cx*axis_x+cy*axis_y+midpoint;
end

end

function tricent=findtrianglecentroid(points)
% centroid of a triangular plane is just the midpoint of vertices. This
% allows faster calculations for triangles.
tricent=sum(points,1)/size(points,1);
end