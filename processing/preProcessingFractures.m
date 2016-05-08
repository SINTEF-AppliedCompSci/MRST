function [G,fracplanes] = preProcessingFractures(G, GlobTri, fracplanes, varargin)

opt = struct('inPolygonTolerance',       0, ...
             'uniqueTolerance'   ,    1e-2, ...
             'pointDensity'      ,    10   );
         
opt = merge_options(opt, varargin{:});
pdens = opt.pointDensity;
%
dims = max(G.nodes.coords);
fracScaled = fracplanes; 
Ar = zeros(numel(fracScaled),1);
Ar_scaled = zeros(numel(fracScaled),1);
for i = 1:numel(fracScaled)
    fracScaled(i).points = [fracScaled(i).points(:,1)./dims(1), ...
                            fracScaled(i).points(:,2)./dims(2), ...
                            fracScaled(i).points(:,3)./dims(3)];
    Ar(i) = polyArea3D(fracplanes(i).points);
    Ar_scaled(i) = polyArea3D(fracScaled(i).points);
end
lx = dims(1); ly = dims(2); lz = dims(3);
scale = 500*min(G.cells.volumes./lx./ly./lz);
esize_scaled = scale/max(Ar_scaled);

% Plane Normals
fracScaled = getPlaneNormals(fracScaled);
fracplanes = getPlaneNormals(fracplanes);
Sign = establishSign(G, fracplanes);

%
[fraCells,remove] = markcells(G, fracplanes, 'Sign', Sign); %, 'GlobTri', GlobTri);
if ~isempty(remove)
    fracplanes = fracplanes(setdiff(1:numel(fracplanes),remove));
    fracScaled = fracScaled(setdiff(1:numel(fracScaled),remove));
end
rectangular = findRectangularFractures(fracplanes);
%
utol = opt.uniqueTolerance;
Gf = struct;

CI = cell(numel(fracplanes),1);
cstart = G.cells.num + 1;
fstart = G.faces.num + 1;
nstart = G.nodes.num + 1;
G.nnc = struct;
G.nnc.cells = [];
G.nnc.CI = [];
G.nnc.area = [];
for i = 1:numel(fracplanes)
    nc = G.nodes.coords(unique(gridCellNodes(G,fraCells{i,1})),:);
    cc = G.cells.centroids(fraCells{i,1},:);
    %
    tripoints = generatePointsOnPlane(fracplanes(i).points, 'normal', ...
        fracplanes(i).normal, 'ratio', pdens);
    tripoints = uniquetol([cc;nc;tripoints], utol, 'ByRows', true);
    %
    tri = delaunayTriangulation(tripoints);
    
    fieldname = ['Frac',num2str(i)];
    isrectangular = ismember(i,rectangular);
    scale = max(dims(3)/dims(1),dims(1)/dims(3));
%     Gf.(fieldname) = gridFractureDistmesh(G, fracplanes(i), fracScaled(i), esize_scaled, ...
%                      'type', 'pebi','rectangular', isrectangular, 'scale', scale);
    Gf.(fieldname) = gridPlanarFracture(G, fracplanes(i), fracScaled(i), esize_scaled*10, ...
                     'type', 1,'rectangular', isrectangular, 'scale', scale);
    Gf.(fieldname).cells.start = cstart;
    Gf.(fieldname).faces.start = fstart;
    Gf.(fieldname).nodes.start = nstart;
    
    CI{i,1} = totalCI(tri,fracplanes(i));
    G = fracMatrixConnections(G, Gf.(fieldname), CI{i,1}, fraCells{i,1}, ...
        polyArea3D(fracplanes(i).points));
    
    cstart = cstart + Gf.(fieldname).cells.num;
    fstart = fstart + Gf.(fieldname).faces.num;
    nstart = nstart + Gf.(fieldname).nodes.num;
end
G.FracGrid = Gf;
return

function ci = totalCI(tri,fracplane)

T = tri.ConnectivityList;
P = tri.Points;

cent = [mean([P(T(:,1),1),P(T(:,2),1),P(T(:,3),1),P(T(:,4),1)],2),...
        mean([P(T(:,1),2),P(T(:,2),2),P(T(:,3),2),P(T(:,4),2)],2),...
        mean([P(T(:,1),3),P(T(:,2),3),P(T(:,3),3),P(T(:,4),3)],2)];

points = [fracplane.points;fracplane.points(1,:)];
sign = zeros(size(cent,1),size(fracplane.boundNormals,1));
tol = eps*100;
for j = 1:size(fracplane.boundNormals,1)
        normal = fracplane.boundNormals(j,:);
        A = normal(1); B = normal(2); C = normal(3);
        D = -dot(normal,mean(points(j:j+1,:),1));
        
        sn = A*cent(:,1) + B*cent(:,2) + C*cent(:,3) + D;
        sn(abs(sn)>tol) = sn(abs(sn)>tol)./abs(sn(abs(sn)>tol));
        sn(isnan(sn)) = 0;
        sn(abs(sn)<tol) = 0;
        sign(:,j) = sn;
end

triToUse = find(sum(sign>=0,2)==size(fracplane.points,1));

P21 = P(T(:,2),:)-P(T(:,1),:);
P31 = P(T(:,3),:)-P(T(:,1),:);
P41 = P(T(:,4),:)-P(T(:,1),:);
V = 1/6*abs(dot(P21,cross(P31,P41,2),2));

normal = fracplane.normal./norm(fracplane.normal);
A = normal(1); B = normal(2); C = normal(3);
D = -dot(normal,points(1,:));

dist = abs((A*cent(triToUse,1) + B*cent(triToUse,2) + C*cent(triToUse,3) + D))./norm(normal);
dV = dist.*V(triToUse);
davg = sum(dV)/sum(V(triToUse));
ci = polyArea3D(fracplane.points)/davg;
return