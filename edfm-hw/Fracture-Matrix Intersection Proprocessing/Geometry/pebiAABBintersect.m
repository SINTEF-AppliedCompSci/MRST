function [truth,area,type,location,points]=pebiAABBintersect(PEBInodes,AABBnodes,tol,varargin)
% Given 3D AABB cell nodes and 2D PEBI nodes, determine if they intersect.
% If they intersect, the area and type of intersections are provided. If no
% intersection, area is 0 and type is 'none'. The output type will only
% take 'polygon' or 'none'. Line and point intersections are regarded as
% type 'none' and truth will be false.

opt=struct('boundaryoneside',false,'planenormal',[0 0 0]);
opt=merge_options(opt,varargin{:});

% Arrange PEBI nodes
normal=getnormal(PEBInodes,tol);
numvertices=size(PEBInodes,1);
PEBInodes=arrangenodes(PEBInodes,normal,'numvertices',numvertices);

% Set up triangles
numtriangles=numvertices-2;
triangles=[ones(numtriangles,1),(2:(numvertices-1))',(3:numvertices)'];

% Check triangle AABB intersection
points_cell=cell(numtriangles,1);
area=0; location=0;
for i=1:numtriangles
    trinodes=PEBInodes(triangles(i,:),:);
    [xtruth,xpoints,xtype,xlocation] = triangleAABBintersect(trinodes,AABBnodes,tol);
    points_cell{i}=xpoints;
    if strcmp(xtype,'polygon') && xtruth
        area=area+convexpolygonarea(xpoints,tol,'normal',normal);
        
        switch xlocation
            case 'none'
                xlocation=0;
            case 'boundary';
                xlocation=1;
            case 'interior'
                xlocation=2;
        end
        location=max(location,xlocation);
    end
end

% Check summed area. If zero, then no intersection, else, intersected.
if area<tol
    truth=false; area=0; type='none'; location='none'; points=[];
    return;
else
    truth=true; type='polygon';
    switch location
        case 1
            location='boundary';
        case 2
            location='interior';
    end
end

points=vertcat(points_cell{:});
points=uniquetol(points,tol,'ByRows',true,'DataScale',1);
points=arrangenodes(points,normal);
points=removeredundantnodes(points,tol);

% pick only cells on positive side of planenormal
if opt.boundaryoneside && strcmp(location,'boundary')
    AABBcentroid=sum(AABBnodes)/size(AABBnodes,1);
    dir=dot(AABBcentroid-PEBInodes(1,:),opt.planenormal);
    dir=dir/abs(dir);
    truth=dir>0;  
end

end

function points=removeredundantnodes(points,tol)
rotpts=[points(end,:);points;points(1,:)];

numpoints=size(points,1);
keep=true(numpoints,1);
for i=1:numpoints
    p0=rotpts(i,:);
    p1=rotpts(i+1,:);
    p2=rotpts(i+2,:);
    
    normal=cross(p1-p0,p1-p2);
    
    if abs(norm(normal))<tol
        keep(i)=false;
    end    
end

points=points(keep,:);

end