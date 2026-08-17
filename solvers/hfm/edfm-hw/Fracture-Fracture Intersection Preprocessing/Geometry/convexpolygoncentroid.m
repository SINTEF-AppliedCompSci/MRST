function centroid=convexpolygoncentroid(nodes,tol)
% Given vertices in 'nodes', find the centroid of the polygon. This only
% works for convex polygons.
%
% REVISED: the original built a non-orthonormal in-plane frame
%   axis_y = cross(cross(axis_x,axis_y_int), axis_x)
% whose axis_y had length sin(theta) != 1, making the centroid depend on the
% input vertex order. Here the in-plane frame is made orthonormal using the
% plane normal, so the centroid is order-independent.

numverts=size(nodes,1);

if numverts<=3
    centroid=findtrianglecentroid(nodes);
    return;
else
    normal=getnormal(nodes,tol);
    nodes=arrangenodes(nodes,normal,'numvertices',numverts);

    % determine an ORTHONORMAL in-plane x-y frame of points
    midpoint=sum(nodes,1)/numverts;
    axis_x=(nodes(1,:)-nodes(3,:)); axis_x=axis_x/norm(axis_x);
    axis_y=cross(normal,axis_x);    axis_y=axis_y/norm(axis_y);

    % in-plane x-y coordinates for the 3D points
    inplanepointsx=zeros(numverts,1);inplanepointsy=zeros(numverts,1);
    for i=1:numverts
        inplanepointsx(i)=dot(nodes(i,:)-midpoint,axis_x);
        inplanepointsy(i)=dot(nodes(i,:)-midpoint,axis_y);
    end

    x=[inplanepointsx;inplanepointsx(1)];
    y=[inplanepointsy;inplanepointsy(1)];

    %generate in-plane centroid location
    cx=0;cy=0;area=0;
    for i=1:numverts
        mult=x(i)*y(i+1)-x(i+1)*y(i);
        cx=cx+(x(i)+x(i+1))*mult;
        cy=cy+(y(i)+y(i+1))*mult;
        area=area+0.5*mult;
    end
    cx=cx/(6*area);
    cy=cy/(6*area);

    % convert in plane coordinates to original coordinates
    centroid=cx*axis_x+cy*axis_y+midpoint;
end
end

function tricent=findtrianglecentroid(points)
% centroid of a triangular plane is just the midpoint of vertices.
tricent=sum(points,1)/size(points,1);
end
