function [truth,line,length,df]=PEBIPEBIintersection3D(P1,P2,tol)
% this function takes in 2D PEBI nodes and checks for intersection between
% the two cells. The points for the line, length of the line and df value
% are returned. The df value is for the purpose of transmissibility
% calculation in fracture-fracture intersections.

% arrange P1 nodes
normalP1=getnormal(P1,tol);
numverticesP1=size(P1,1);
P1=arrangenodes(P1,normalP1,'numvertices',numverticesP1);

% arrange P2 nodes
normalP2=getnormal(P2,tol);
numverticesP2=size(P2,1);
P2=arrangenodes(P2,normalP2,'numvertices',numverticesP2);

% set up P1 triangles
numtrianglesP1=numverticesP1-2;
trianglesP1=[ones(numtrianglesP1,1),(2:(numverticesP1-1))',(3:numverticesP1)'];

% set up P2 triangles
numtrianglesP2=numverticesP2-2;
trianglesP2=[ones(numtrianglesP2,1),(2:(numverticesP2-1))',(3:numverticesP2)'];

% test all triangle combinations
line=[];
for i=1:numtrianglesP1
    P1Ti=P1(trianglesP1(i,:),:); % nodes for i-th triangle for P1
    for j=1:numtrianglesP2
        P2Tj=P2(trianglesP2(j,:),:); % nodes for the j-th triangle for P2
        [xtruth,xline,xlength,~]=triangletriangleintersection3D(P1Ti,P2Tj,tol);
        if xlength>tol && xtruth
            line=collinearlineunion(line,xline);
        end
    end
end

% if line is an empty set (meaning none of the triangles intersect), there
% is no intersection line.
if isempty(line)
    truth=false; length=0; df=[];
else
    truth=true;
    length=norm(line(1,:)-line(2,:));
    df=calcPEBIdf(P1,P2,line,tol);
end

end