function F=fractureplanegeneralgrid(V,C,planepoints,planenormal,aperture,tol)
% Generates a Fracture grid following grid_structure, given V and C where V
% is a list of vertex coordinates (3D and repeating allowed) and C is a
% cell for which each element represents a cell and contains the indices of
% V which correspond to its vertices. The fracture plane vertices and
% normal direction are also required.

% V = round(V*10^9)/10^9;

%% Transformation to a 2D set of points
% Rotate points about the origin to align with the xy plane
unitNormal = planenormal/norm(planenormal);
xyp_plane = rotatePlane(planepoints,[0 0 1],'useNormal',unitNormal);
xyp = rotatePlane(V,[0 0 1],'useNormal',unitNormal);

% Extract xy-coordinates and store the z coordinate
xyp_plane=xyp_plane(:,1:2);

z=xyp(1,3);
xyp = xyp(:,1:2);

% Rotate axis to align 1 edge of the planar polygon with the x-axis 
diffp = diff(xyp_plane);
theta = -atan(diffp(1,2)/diffp(1,1));
R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
xyp = transpose(R*xyp');
V=xyp;



%% Create grid
% [V,~,IC] = unique(V,'rows');
[V,~,IC]=uniquetol(V,tol,'ByRows',true,'DataScale',1);

% figure;
% scatter(V(:,1),V(:,2));

F.cells.num = numel(C);
F.cells.facePos = cumsum([1; cellfun(@numel, C)]);
C2 = cellfun(@(c) [c,c(1)],C,'un', false);
faces = horzcat(C2{:})';
faces = [faces(1:end-1), faces(2:end)];
addTo = cumsum([0;ones(F.cells.num-2,1)]);

faces(F.cells.facePos(2:end-1) + addTo,:) = [];
faces = sort(IC(faces),2);

[nodes,~,faces] = unique(faces, 'rows');

F.cells.faces = faces;

F.nodes.coords = V;
F.nodes.num = size(V,1);

F.faces.num = max(faces);
F.faces.nodes = reshape(nodes',[],1);
F.faces.nodePos = (1:2:2*F.faces.num+1)';


cellNo            = rldecode(1:F.cells.num, diff(F.cells.facePos), 2).';
F.faces.neighbors = zeros(F.faces.num,2);
for i = 1:F.faces.num
    neigh = F.cells.faces==i;
    if sum(neigh)==2
        F.faces.neighbors(i,[1,2]) = cellNo(neigh);
    else
          F.faces.neighbors(i,:) = [cellNo(neigh),0];
    end
end

F.griddim = 2;
F.type = {'blah'};
F = sortEdges(F);

%% Transform back to original orientation and make 3D grid
% Rotate axis back to original position
nc = F.nodes.coords;
R = [cos(theta) sin(theta); -sin(theta) cos(theta)];
nc = transpose(R*nc');
F.nodes.coords = nc;

% Make Layered Grid
F = computeGeometry(makeLayeredGrid(F,1));

% Add z-coordinate of the original set of points post rotation
nc = [nc,repmat(z,size(nc,1),1)];

% Rotate points about the origin back to their original position
use_points = rotatePlane(nc,unitNormal,'useNormal',[0 0 1]);

% Expand along normal to create an aperture
tdist = aperture/2;
pminus = use_points - repmat(unitNormal,size(use_points,1),1)*tdist; %correct
pplus = use_points + repmat(unitNormal,size(use_points,1),1)*tdist; %correct
nc = [pminus;pplus];

% Add new node coordinates to fracture grid and re-compute geometry
F.nodes.coords = nc;
F = computeGeometry(F);
F.cells.volumes = abs(F.cells.volumes);
F.faces.areas = abs(F.faces.areas);
F.type='General';
F.nodes2D.coords = use_points;


end