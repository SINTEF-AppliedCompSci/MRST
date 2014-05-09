function G = computeGeometryVE(G, varargin)
%Compute geometry of grid of top-surface grid
%
% SYNOPSIS:
%   G = computeGeometryVE(G)
%   G = computeGeometryVE(G, 'pn', pv, ...)
%
% PARAMETERS:
%   G       - Grid structure as defined by function 'topSurfaceGrid'.
%
%   'pn'/pv - List of 'key'/value pairs for supplying optional parameters.
%             The supported options are
%               - Verbose -- Whether or not to display verbose output as
%                       the process progresses.  Possible values are
%                       TRUE and FALSE.  Default value equals mrstVerbose.
%
% RETURNS:
%   G - Grid structure with added fields:
%         - cells
%             - volumes
%             - centroids
%             - z
%         - faces
%             - areas
%             - normals                !!! AREA WEIGHTED !!!
%             - centroids
%             - z

%
% COMMENTS:
%   PLEASE NOTE: Face normals have length equal to face areas

%{
#COPYRIGHT#
%}

% $Date: 2012-01-30 11:41:03 +0100 (Mon, 30 Jan 2012) $
% $Revision: 9020 $

% Setup
assert(size(G.faces.nodes, 2)==1);
opt     = struct('verbose', mrstVerbose);
opt     = merge_options(opt, varargin{:});

if size(G.nodes.coords,2)==2
   coords    = [G.nodes.coords G.nodes.z];
else
   assert(false);
end

dispif(opt.verbose, 'Computing normals, areas, and centroids...\t');
ticif (opt.verbose)

faceEdges = reshape(G.faces.nodes,2,[])';

edgeLength(:,:) = coords(faceEdges(:,2),:) - coords(faceEdges(:,1),:);
faceAreas     = sqrt(sum(edgeLength.*edgeLength,2));
faceCentroids = (coords(faceEdges(:,2),:)+ coords(faceEdges(:,1),:))/2;
faceNormals   = [edgeLength(:,2),-edgeLength(:,1)];

tocif(opt.verbose)
dispif(opt.verbose, 'Computing cell volumes and centroids ...\t\t');
ticif (opt.verbose)

cellno = rldecode(1:G.cells.num, diff(G.cells.facePos), 2)';

numfaces=diff(G.cells.facePos);

cCenter       = zeros(G.cells.num, 3);
cCenter(:,1)  = accumarray(cellno, faceCentroids(G.cells.faces(:,1), 1));
cCenter(:,2)  = accumarray(cellno, faceCentroids(G.cells.faces(:,1), 2));
cCenter(:,3)  = accumarray(cellno, faceCentroids(G.cells.faces(:,1), 3));
cCenter       = bsxfun(@rdivide, cCenter, double(numfaces));

a = coords(faceEdges(G.cells.faces(:,1),2),:) - cCenter(cellno,:);
b = coords(faceEdges(G.cells.faces(:,1),1),:) - cCenter(cellno,:);
subArea      = sqrt(sum(cross(a,b,2).^2,2))*0.5;
subCentroid  = bsxfun(@plus, cCenter(cellno,:), 2*faceCentroids(G.cells.faces(:,1),:))/3;
cellVolumes  = accumarray(cellno, subArea);

cellCentroids      = zeros(G.cells.num, 3);
cellCentroids(:,1) = accumarray(cellno, bsxfun(@times, subArea, subCentroid(:,1)));
cellCentroids(:,2) = accumarray(cellno, bsxfun(@times, subArea, subCentroid(:,2)));
cellCentroids(:,3) = accumarray(cellno, bsxfun(@times, subArea, subCentroid(:,3)));
cellCentroids      = bsxfun(@rdivide, cellCentroids, cellVolumes);
tocif(opt.verbose)

% Update grid
G.faces.areas     = faceAreas;
G.faces.normals   = faceNormals;
G.faces.centroids = faceCentroids(:,1:2);
G.faces.z = faceCentroids(:,3);

G.cells.volumes   = cellVolumes;
G.cells.centroids = cellCentroids(:,1:2);
G.cells.z = cellCentroids(:,3);
G.type = [G.type, mfilename];

