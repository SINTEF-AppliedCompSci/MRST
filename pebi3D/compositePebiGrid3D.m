function [G,F] = compositePebiGrid3D(celldim, pdim, varargin)
% Construct a 3D composite Pebi grid. This function creates a cartesian 
% background grid with pebi-cells around face constraints and cell constraints.
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
%   faceConstraints   - OPTIONAL.
%                       Default value empty. A nfx1 cell array of 
%                       triangulations. Each triangulation should be a 2D
%                       surface embedded in 3D. The faces of the grid
%                       will conform to the surfaces.
%   FCRho         - OPTIONAL.
%                       Default value sqrt(sum(pdims./celldim.^2)). 
%                       FCRho is a cell array of length nf or 1.
%                       Each element in FCRho is a function. 
%                       The function must map from a nx3 array to nx1 
%                       array. FCRho{i}(faceConstraints{i}.Points) define the 
%                       radius of the balls placed along surface i. If 
%                       FCRho has one element, this function is applied
%                       to all surfaces.
%   cellConstraints        - OPTIONAL.
%                       Default value {{}}. A nwx1 cell array of vectors. 
%                       Each element, of size kx3, contains the coordinates 
%                       of a line. The lines is assumed to be linear
%                       between the coordinates. A set of sites will be
%                       placed along these lines.
%   CCRho           - OPTIONAL.
%                       Default value sqrt(sum(pdims./celldim.^2)). CCRho
%                       is a cell array of length nwx1 or 1. Each element
%                       is a function. The fucntion must map from a nx3
%                       array to nx1 array. CCRho{i}(p) defines the
%                       desired distance between the sites that are placed
%                       along the lines given by cellConstraints{i}. If 
%                       CCRho has one element, this function is applied to
%                       all well paths.
% RETURNS:
%   G                - Valid grid definition.  
%                        Contains extra fields:
%                          - G.cells.tag is TRUE for all cell constraint cells.
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
if ~all(celldim > 0)
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
             'mlqtMaxLevel',    0,    ...
             'mlqtLevelSteps',  -1,   ...
             'faceConstraints', {{}},      ...
             'FCRho',  {{@(x) sqrt(pdim./sum(celldim.^2))*ones(size(x,1),1)}});          
opt = merge_options(opt, varargin{:});

FCRho = opt.FCRho;
mlqtMaxLevel   = opt.mlqtMaxLevel;
mlqtLevelSteps = opt.mlqtLevelSteps;

if numel(FCRho) == 1
  FCRho = repmat(FCRho,numel(opt.faceConstraints),1);num2cell(FCRho, 2);
else
  assert(numel(FCRho) == numel(opt.faceConstraint),...
    'Number of FCRho must either be 1 or numel(faceConstraint)');
end
CCRho = opt.CCRho;
if numel(CCRho) == 1
  CCRho = repmat(CCRho,numel(opt.cellConstraints),1);num2cell(CCRho, 2);
else
  assert(numel(CCRho) == numel(opt.cellConstraints),...
    'Number of CCRho must either be 1 or numel(cellConstraints)');
end

if numel(FCRho) == 1
  FCRho = repmat(FCRho,numel(opt.faceConstraints),1);num2cell(FCRho, 2);
else
  assert(numel(FCRho) == numel(opt.faceConstraints),...
    'Number of FCRho must either be 1 or numel(faceConstraint)');
end

assert(mlqtMaxLevel>=0, 'mlqtMaxLevel must be greater or equal 0');

% Create face constraint sites
F = surfaceSites3D(opt.faceConstraints,FCRho);

% Remove conflict points at fault intersections
F = removeSurfaceConflictSites3D(F);

% Create cell constraint sites
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

% Possible refine mesh
if ~isempty(W.pts)
    varArg = {'level', 1, 'maxLev', mlqtMaxLevel, 'distTol', mlqtLevelSteps};
    res = {};
    for i = 1:size(rSites, 1)
        res = [res; mlqt(rSites(i,:), W.pts , [dx, dy , dz], varArg{:})];
    end
    rSites = vertcat(res{:, 1});
end
% Remove conflict sites
rSites = lineSufCond3D(rSites, W);
rSites = surfaceSufCond3D(rSites,F.c.CC,F.c.R);
W.pts  = surfaceSufCond3D(W.pts,F.c.CC,F.c.R);

% Create grid
pts = [F.f.pts; W.pts; rSites];
%bdrDT = delaunayTriangulation(bdr);
%G = clippedPebi3D(pts,bdrDT);
G = mirroredPebi3D(pts, bdr);

% Tag well cells 
G.cells.tag = false(G.cells.num,1);
G.cells.tag(size(F.f.pts,1)+1:size(F.f.pts,1)+size(W.pts,1)) = true;

end
