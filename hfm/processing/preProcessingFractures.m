function [G,fracplanes] = preProcessingFractures(G, fracplanes, varargin)
% preProcessingFractures identifies fracture-matrix connections, imposes a
% fracture grid and computes the fracture-matrix conductivity index.
%
% SYNOPSIS:
%   [G,fracplanes] = preProcessingFractures(G, fracplanes)
%   [G,fracplanes] = preProcessingFractures(G, fracplanes, 'pn1', pv1)
%
% REQUIRED PARAMETERS:
%
%   G           - Grid data structure containing geometrical information.
%
%   fracplanes  - 1-by-n structure where n is the number of fractures. Each
%                 column contains a set of coplanar points that define
%                 the fracture polygon and a value for the fracture
%                 aperture.
%
% OPTIONAL PARAMETERS:
%   
%   fractureGridType   - Type of fracture mesh desired. Possible options:
%                        (a) 1 - Cartesian mesh: Possible only when the
%                        fracture plane is rectangular.
%
%                        (b) 2 - Triangle mesh: Uses the external library
%                        'distmesh' to create the nodes and connectivity
%                        list required to create a triangular grid.
%
%                        (c) 3 - PEBI mesh: Voronoi grid. External
%                        library 'distmesh' required to create underlying
%                        triangulation.
%
%   fractureCellSize   - Dimensionless element size (>0 and <1) for the
%                        fracture grid.
%                
%   GlobTri            - See globalTriangulation. 
%
%   pointDensity       - A measure to control the number density of points
%                        generated on a fracture plane for processing
%                        purposes.
%
%   inPolygonTolerance - Tolerance for checking if a set of points lie
%                        inside a polygon.
%
%   uniqueTolerance    - Tolerance for the matlab function unique.
%
% RETURNS:
%   G           - Grid data structure with added fields G.nnc and
%                 G.FracGrid. G.FracGrid contains the grid information for
%                 each fracture. G.nnc contains relevant information
%                 about fracture matrix connections, which are defined as
%                 NNC's to compute transmissibilities.
%
%   fracplanes  - Structure with added fields normal and boundNormals for
%                 each fracture polygon.
%
% SEE ALSO:
%   markcells, getPlaneNormals, gridPlanarFracture, fracMatrixConnections


%{
Copyright 2009-2015: TU Delft and SINTEF ICT, Applied Mathematics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}

opt = struct('fractureGridType'  ,    1     , ...
             'fractureCellSize'  ,    0     , ...
             'pointDensity'      ,    10    , ...
             'inPolygonTolerance',    0     , ...
             'uniqueTolerance'   ,    1e-2  , ...
             'GlobTri'           ,    struct('Tri', [], 'map', []));
         
opt = merge_options(opt, varargin{:});
if ~any(opt.fractureGridType == [1,2,3])
    fprintf('\nWrong input for fracture grid type, selecting cartesian grid type.\n');
    opt.fractureGridType = 1;
end

pdens = opt.pointDensity;
%
dims = max(G.nodes.coords);

% dimensionless fracture plane definition for gridding purposes
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

% Define element size for gridding
if opt.fractureCellSize == 0
    esize_scaled = scale/max(Ar_scaled);
else
    esize_scaled = opt.fractureCellSize;
end

% Plane normals
fracScaled = getPlaneNormals(fracScaled);
fracplanes = getPlaneNormals(fracplanes);
Sign = establishSign(G, fracplanes);

% Mark cells containing fractures
[fraCells,remove] = markcells(G, fracplanes, ...
                    'Sign', Sign, ...
                    'GlobTri', opt.GlobTri, ...
                    'inPolygonTolerance', opt.inPolygonTolerance);
if ~isempty(remove)
    fracplanes = fracplanes(setdiff(1:numel(fracplanes),remove));
    fracScaled = fracScaled(setdiff(1:numel(fracScaled),remove));
end
rectangular = findRectangularFractures(fracplanes);
%
utol = opt.uniqueTolerance;
Gf = struct;
dispif(mrstVerbose, 'Computing CI and constructing fracture grid...\n\n');
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
    tri = delaunayTriangulation(tripoints);
    %
    fieldname = ['Frac',num2str(i)];
    isrectangular = ismember(i,rectangular);
    scale = max(dims(3)/dims(1),dims(1)/dims(3));
    Gf.(fieldname) = gridPlanarFracture(G, fracplanes(i), fracScaled(i), ...
                     'rectangular', isrectangular, 'scale', scale, ...
                     'gridType', opt.fractureGridType, ...
                     'cellSize', esize_scaled);
    Gf.(fieldname).cells.start = cstart;
    Gf.(fieldname).faces.start = fstart;
    Gf.(fieldname).nodes.start = nstart;
    %
    CI{i,1} = totalCI(tri,fracplanes(i));
    G = fracMatrixConnections(G, Gf.(fieldname), CI{i,1}, fraCells{i,1}, ...
        polyArea3D(fracplanes(i).points), 'GlobTri', opt.GlobTri);
    
    cstart = cstart + Gf.(fieldname).cells.num;
    fstart = fstart + Gf.(fieldname).faces.num;
    nstart = nstart + Gf.(fieldname).nodes.num;
end
G.FracGrid = Gf;
return

function ci = totalCI(tri,fracplane)
% Computes the conductivity index for the entire fracture surface given the
% matrix cells it is connected with.

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