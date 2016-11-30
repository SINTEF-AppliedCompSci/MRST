function [G,dt,symV] = clippedPebi3D(p, bound)
% Create a clipped Pebi Grid
%
% SYNOPSIS:
%   G = clippedPebi3D(p,bound)
%
% PARAMETERS:
%   p         A nx3 array containing the voronoi sites. 
%   bound     A Delaunay triangulation of the domain. The returned grid
%             will be restricted to this domain. The triangulation must
%             contain the following elements, as described in the matlab
%             routine delaunayTriangulation:
%               - bound.Points
%               - bound.freeBoundary
%               - bound. pointLocation
%
% RETURNS:
%   G                - Valid MRST grid definition.  
%
% EXAMPLE:
% dt = 0.15;
% [X,Y,Z] = ndgrid(dt:dt:1-dt);
% X(1:2:end) = X(1:2:end) + dt/2;
% Y(1:2:end) = Y(1:2:end) + dt/2;
% Z(1:2:end) = Z(1:2:end) + dt/2;
% pts = [X(:), Y(:),Z(:)];
% boundPts = [0,0,0; 1,0,0; 1,1,0; 0,1,0; 0,0,1; 1,0,1; 1,1,1;0,1,1];
% bound = delaunayTriangulation(boundPts);
% G = clippedPebi3D(pts, bound);
% plotGrid(G);
% axis equal
%
% SEE ALSO:
%   compositePebiGrid, pebi, voronoi2mrst, clipGrid.

%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2016 Runar Lie Berge. See COPYRIGHT.TXT for details.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}

assert(size(p, 2) == 3, ...
      ['Function ''%s'' is only supported in three ', ...
       'space dimensions.'], mfilename);

dt = delaunayTriangulation(p);

dtB.ConnectivityList = bound.freeBoundary;
dtB.Points = bound.Points;

% Clip grid against boundary
[Vext,Cext,symV] = clipGrid(dt,dtB);


% Add inner vertices, i.e. intersection of three bisections
[VInt, CInt] = voronoin(p);%dt.voronoiDiagram();

% Remove points outside domain
keep = [false;~isnan(bound.pointLocation(VInt(2:end,:)))];

VInt = VInt(keep,:);
keepNum = find(keep);
newIdx  = zeros(size(VInt,1),1);
newIdx(keepNum) = (1:size(VInt,1))';
CInt = cellfun(@(c) newIdx(intersect(c,keepNum))', CInt,'uniformOutput',false);
ni   = size(VInt,1);
V    = [VInt;Vext];
symV = [cell(size(VInt,1),1);symV];
C    = cellfun(@(cint,cext) [cint,cext+ni], CInt, Cext,'uniformOutput',false);

[V,IA,IC] = uniquetol(V,100*eps,'byRows',true);
symV = symV(IA);
C    = cellfun(@(c) unique(IC(c))', C,'UniformOutput',false);

% Tansfrom to Mrst grid
G    = voronoi2mrst(V,C);
end