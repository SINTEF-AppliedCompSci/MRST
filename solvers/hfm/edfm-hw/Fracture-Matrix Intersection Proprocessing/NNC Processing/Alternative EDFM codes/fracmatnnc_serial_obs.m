function [cells,CI,area]=fracmatnnc_serial(F,Gm,tol)
% INDIVIDUAL FRACTURE-MATRIX NNC GENERATOR - This function takes in the
% fracture grid of an individual fracture and the matrix grid. It then
% generates outputs cells (which contains global indices in each row for
% every frac-mat nnc), CI (which contains the CI for each frac-mat nnc,
% corresponding to the cells array) and area (which contains the area for
% each intersection, corresponding to the cells array).

% METHODOLOGY - We first get the plane normal for the fracture and do a
% quick filter using funciton checkplaneAABBintersect. Then, we save the
% indices of the cells that could possibly intersect the fracture.
%
% Then, for every fracture gridcell (which is triangular), we test for
% triangle-AABB intersection with every possible intersecting matrix cell
% using the triangleAABBintersect function.
%
% If an intersection is detected, and the intersection is a polygon (which
% has an area), we go on to calculate polygon area and average normal
% distance (per Moinfar et al, 2013). The CI is then calculated, global
% indices calculated and all data appended to cells, CI and area.

%% FILTERING STEP
nummatcells=Gm.cells.num;

% find normal of fracture plane using one of the fracture gridcells
% (triangular)
[nodeind,~] = gridCellNodes(F,1);
nodes=F.nodes.coords(nodeind,:);
dir1=nodes(2,:)-nodes(1,:);
dir2=nodes(3,:)-nodes(1,:);
planenormal=cross(dir1,dir2); planenormal=planenormal/norm(planenormal); % make unit normal
planepoint=nodes(1,:);

possiblemcells=-1*ones(nummatcells,1); % preallocate max number of cells, remove excess later
count=0;
[cn,cpos]=gridCellNodes(Gm,1:nummatcells);
for i=1:nummatcells
    % find out if cell 'i' intersects infinite plane
    mcellnodeind=cn(cpos(i):(cpos(i+1)-1));
    mcellnodes=Gm.nodes.coords(mcellnodeind,:);
    [xtruth,xtype]=checkplaneAABBintersect(mcellnodes,planepoint,planenormal,tol);
    if xtruth && strcmp(xtype,'interior')
        count=count+1;
        possiblemcells(count)=i;
    end    
end
possiblemcells=removeexcess(possiblemcells,-1);
% at this point, we have an array possiblecells which contains indices of
% cells that could potentially intersect the fracture. We proceed to a more
% rigorous intersection detection in the next step.

%% RIGOROUS STEP
numpossiblemcells=length(possiblemcells);

fcellstart=F.cells.start;
numfraccells=F.cells.num;

[fn,fpos]=gridCellNodes(F,1:numfraccells);

maxallocate=numfraccells*numpossiblemcells;
cells=-1*ones(maxallocate,2);CI=-1*ones(maxallocate,1);area=-1*ones(maxallocate,1); % preallocate max, remove excess later
count=0;
davg=ones(numpossiblemcells,1)*(-1);
for i=1:numfraccells
    disp(['Processing fracture cell [',num2str(i),'/',num2str(numfraccells),']...']);
    globalfcellind=i+fcellstart-1;
    fcellnodeind=fn(fpos(i):(fpos(i+1)-1));
    fcellnodes=F.nodes.coords(fcellnodeind,:); % 6x3 matrix, each row containing xyz coords of nodes
    fcellnodes=0.5*(fcellnodes(1:3,:)+fcellnodes(4:6,:));
    if degenerate(fcellnodes,tol) % some triangles have 0 area. Meshing algorithm not strong enough.
        continue;
    end
    
    for j=1:numpossiblemcells
        cellind=possiblemcells(j);
        mcellnodeind=cn(cpos(cellind):(cpos(cellind+1)-1));
        mcellnodes=Gm.nodes.coords(mcellnodeind,:); % 8x3 matrix, each row containing xyz coords of nodes
        [xtruth,xpoints,xtype,~]=triangleAABBintersect(fcellnodes,mcellnodes,tol);
        
        if xtruth && strcmp(xtype,'polygon')
            % if fracture cell intersects a matrix cell, append details to
            % cells, area and CI.
            count=count+1;
            cells(count,:)=[cellind,globalfcellind];
            
            xarea=convexpolygonarea(xpoints,tol,'normal',planenormal); % implicit assumption here is that intersections are all convex
            area(count)=xarea;
            [xCI,davg(j)]=calcfracmatCI(mcellnodes,xarea,planenormal,planepoint,davg(j),tol);
            CI(count)=xCI;
        end            
    end
end

% once the function gets out of this for block, every fracture grid cell
% would have been compared with every matrix grid cell for intersection. If
% intersection was detected, the global indices, area and CI were saved.

cells=removeexcess(cells,-1);
CI=removeexcess(CI,-1);
area=removeexcess(area,-1);

end

function truth=degenerate(nodes,tol)
% check if a triangle is degenerate. Use heron's formula to check area of
% triangle. If zero, then degenerate.

area=heronsformula(nodes(1,:),nodes(2,:),nodes(3,:),tol);

if abs(area)<tol
    truth=true;
else
    truth=false;
end

end

