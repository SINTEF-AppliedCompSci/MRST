function [G,F] = compositePebiGrid3D(celldim, pdim, varargin)
% Construct a 3D composite Pebi grid. This function creates a cartesian 
% background grid with pebi-cells around faults and wells.
%
% SYNOPSIS:
%   G = compositePebiGrid3D(celldim, pdims)
%   G = compositePebiGrid3D(...,'Name1',Value1,'Name2',Value2,...)
%
% PARAMETERS:
%   celldim           - [xSize,ySize, zSize] the number of cells in each 
%                       coordinate direction.
%   pdim              - Vector, length 2, [xmax, ymax, zmax], of physical
%                       size in units of meters of the computational domain 
%
%   faultSurf         - OPTIONAL.
%                       Default value empty. A nfx1 cell array of 
%                       triangulations. Each triangulation should be a 
%                       surface embedded in 3D. The grid will conform to 
%                       the surfaces.
%   FCRho         - OPTIONAL.
%                       Default value sqrt(sum(pdims./celldim.^2)). 
%                       FCRho is a cell array of length nf or 1.
%                       Each element in FCRho is a function. 
%                       The function must map from a nx3 array to nx1 
%                       array. FCRho{i}(faultSurf{i}.Points)define the 
%                       radius of the balls placed along fault i. If 
%                       FCRho has one element, this function is applied
%                       to all fault surfaces.
%   cellConstraints        - OPTIONAL.
%                       Default value {{}}. A nwx1 cell array of vectors. 
%                       Each element, of size kx3, contains the coordinates 
%                       of a well-trace. The well is assumed to be linear
%                       between the coordinates.
%   CCRho           - OPTIONAL.
%                       Default value sqrt(sum(pdims./celldim.^2)). WellRho
%                       is a cell array of length nwx1 or 1. Each element
%                       is a function. The fucntion must map from a nx3
%                       array to nx1 array. CCRho{i}(p) defines the
%                       desired distance between well sites that are placed
%                       along well path cellConstraints{i}. If CCRho has one
%                       element, this function is applied to all well
%                       paths.
% RETURNS:
%   G                - Valid grid definition.  
%                        Contains extra fields:
%                          - G.cells.tag is TRUE for all well cells.
%       
% EXAMPLE:
%   
%
% SEE ALSO:
%   compositePebiGrid2D, pebi, surfaceSites2D, lineSites2D.

%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2016 Runar Lie Berge. See COPYRIGHT.TXT for details.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}  

% Check input
if ~all(celldim > 0),
   error('CELLDIM must be positive');
end
if numel(celldim)~=3
  error('CELLDIM must have 3 elements')
end

if ~all(pdim>0)
  error('PDIMS must be positive');
end
if numel(pdim)~=3
  error('PDIMS must have 3 elements')
end

% Set options
opt = struct('cellConstraints', {{}}, ...
             'CCRho',   {{@(x) sqrt(pdim./sum(celldim.^2))*ones(size(x,1),1)}},...
             'faultSurf', {{}},      ...
             'FCRho',  {{@(x) sqrt(pdim./sum(celldim.^2))*ones(size(x,1),1)}});          
opt = merge_options(opt, varargin{:});

FCRho = opt.FCRho;
if numel(FCRho) == 1
  FCRho = repmat(FCRho,numel(opt.faultSurf),1);num2cell(FCRho, 2);
else
  assert(numel(FCRho) == numel(opt.faultSurf),...
    'Number of FCRho must either be 1 or numel(faultSurf)');
end
CCRho = opt.CCRho;
if numel(CCRho) == 1
  CCRho = repmat(CCRho,numel(opt.cellConstraints),1);num2cell(CCRho, 2);
else
  assert(numel(CCRho) == numel(opt.cellConstraints),...
    'Number of CCRho must either be 1 or numel(cellConstraints)');
end


% Create fault sites
F = surfaceSites3D(opt.faultSurf,FCRho);

% Remove conflict points at fault intersections
F = removeSurfaceConflictSites3D(F);

% Create well points
W = lineSites3D(opt.cellConstraints, CCRho);

% Create reservoir sites
x = pdim(1); y = pdim(2); z = pdim(3);

dx = x/celldim(1);
dy = y/celldim(2);
dz = z/celldim(3);

xa = linspace(dx/2, x-dx/2, celldim(1)+1);
ya = linspace(dy/2, y-dy/2, celldim(2)+1);
za = linspace(dz/2, z-dz/2, celldim(3)+1);

bdr   = [ 0, 0, 0;  ...
          x, 0, 0;  ...
          x, y, 0;  ...
          0, y, 0;  ...
          0, 0, z;  ...
          x, 0, z;  ...
          x, y, z;  ...
          0, y, z];

[X,Y,Z] = ndgrid(xa,ya,za);
rSites = [X(:), Y(:), Z(:)];

% Remove conflict points
rSites = lineSufCond3D(rSites, W);
rSites = surfaceSufCond3D(rSites,F.c.CC,F.c.R);
W.pts  = surfaceSufCond3D(W.pts,F.c.CC,F.c.R);

% Create grid
pts = [F.f.pts;W.pts;rSites];
%bdrDT = delaunayTriangulation(bdr);
%G = clippedPebi3D(pts,bdrDT);
G = mirroredPebi3D(pts, bdr);

% Tag well cells 
G.cells.tag = false(G.cells.num,1);
G.cells.tag(size(F.f.pts,1)+1:size(F.f.pts,1)+size(W.pts,1)) = true;

end
