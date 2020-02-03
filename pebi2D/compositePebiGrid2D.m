function [G,Pts,F] = compositePebiGrid2D(celldim, pdims, varargin)
% Construct a 2D composite Pebi grid. A cartesian background grid that is 
% refined around face constraints and cell constraints.
%
% SYNOPSIS:
%   G = compositePebiGrid2D(resGridSize, pdims)
%   G = compositePebiGrid2D(...,'Name1',Value1,'Name2',Value2,...)
%
% PARAMETERS:
%   resGridSize       - [xSize,ySize] Size of the reservoir grid cells in x
%                       and y direction.
%
%   pdims             - Vector, length 2, [xmax, ymax], of physical size in
%                       units of meters of the computational domain. 
%
%   cellConstraints - OPTIONAL.
%                     Default value empty. A struct of vectors. Each
%                     vector, size nw x 2, is the coordinates of a
%                     line and the function will place a set of sites along
%                     this line. The lines are assumed to be linear between
%                     the coorinates. If the vector only contains one
%                     coordinate, the line is treated as a point constraint.
%
%   interpolateCC     - OPTIONAL.
%                     Default value is a boolean false value, but the user
%                     can supply an individual value per cell constraint line. 
%                     If false, each segment in the corresponding
%                     constraint line will be represented by at least one
%                     cell. If true, the routine will interpolate along the
%                     curve, which means that sites will not
%                     necessarily fall exactly on the prescribed curve.
%
%   CCFactor        - OPTIONAL.
%                     Default value is 1. This gives the relative grid
%                     size of the cell constraint grid cells compared to
%                     reservoir grid cells. If CCFactor=0.5 the constrained
%                     cells will be about half the size of the reservoir cells.
%
%   mlqtMaxLevel    - OPTIONAL.
%                     Default value 0. Number of refinement steps around 
%                     cell constraints. 
%
%   mlqtLevelSteps  - OPTIONAL.  
%                     Default value -1. Vector of length mlqtMaxLevel 
%                     which specifies the radius of each refinement level.
%                     The default value -1 calls the default level step
%                     in the mlqt function.
%
%   CCRho         - OPTIONAL
%                     Default value @(x) ones(size(x,1),1). Function gives
%                     the relative distance between cell constraint cells. If
%                     CCRho=0.5 in an area the distance between
%                     constrained-cells will be 0.5*CCFactor*resGridSize
%
%   protLayer         - OPTIONAL.
%                       Default set to false. If set to true a protection layer
%                       is added on both sides of the well
%                       .          .             .  Protection Layer
%                                                   protD(dist between points)
%                       .----------.-------------.  Well path
%     
%                       .          .             .  Protection Layer
%     
%   protLayer       - OPTIONAL.
%                     Default set to false. If set to true, a protection layer
%                     is added on both sides of the cell Constraints
%                     .          .             .  Protection Layer
%                                                 protD(dist between sites)
%                     .----------.-------------.  Constraint path
%
%                     .          .             .  Protection Layer
%
%
%   faceConstraints - OPTIONAL
%                     Default value empty. A struct of vectors.  Each
%                     vector, size nf x 2, is the coordinates of a
%                     surface-trace. The surface is assumed to be linear
%                     between the coorinates. The function will place sites
%                     such that the surface is traced by faces of the grid.
%
%   FCFactor        - OPTIONAL.
%                     Default value is 0.5. This gives the relative grid
%                     size of the surface grid cells compared to reservoir
%                     grid cells. If FCFactor=0.5 the surface cells
%                     will be about half the size of the reservoir cells.
%
%   interpolateFC   - OPTIONAL.
%                     Default value is a boolean false value, but the
%                     user can supply an individual value per surface.
%                     If false, each segment in the corresponding surface
%                     curve will be represented by at least one cell. If
%                     true, the routine will interpolate along the
%                     curve, which means that faces will not
%                     necessarily fall exactly on the prescribed curve.
%
%   circleFactor    - OPTIONAL.
%                     Default value 0.6.  Valid values are between 0.5
%                     and 1. circleFactor controll the size of the
%                     circles used to create the surface sites. The
%                     circleFactor is the ratio between the radius
%                     and distance between the circles. A small value will
%                     place the surfaces points close the the surfaces, while
%                     a large value will place the far from the surfaces.
%
%   polyBdr         - OPTIONAL
%                     Default value []. plyBdr is a array of size [k,2].
%                     if k>=3 polyBdr gives the vertices of the reservoir
%                     boundary. For this domain:
%                                .(x1,y1)
%                               / \
%                       (x3,y3).---.(x2,y2)
%                     polyBdr would be [x1,y1;x2,y2;x3,y3]. The set of
%                     vertices must run clockwise or counterclockwise.
%                     Note that if k>=3 the innput pdims will have no
%                     effect.
%
% RETURNS:
%   G              - Valid grid definition.
%                      The fields
%                        - G.cells.tag is TRUE for all cell constraints.
%                        - G.faces.tag is TRUE for all face constraints.
%   Pts            - Array [G.cells.num x 3] of the Pebi sites.
%   F              - Struct with elements as returned from
%                    surfaceSites2D
%
% EXAMPLE:
%   fl = {[0.2,0.2;0.8,0.8]};
%   wl = {[0.2,0.8;0.8,0.2]};
%   G  = compositePebiGrid2D([1/10,1/10],[1,1],'cellConstraints',wl,'faceConstraints',fl)
%   cla, plotGrid(G)
%
% SEE ALSO:
%   compositePebiGrid2D, pebi, surfaceSites2D, lineSites2D.

%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2015-2020 Runar Lie Berge. See COPYRIGHT.TXT for details.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}  

% Set options
opt = struct('cellConstraints', {{}}, ...
             'CCFactor',        1,  ...
             'interpolateCC',   false, ...
             'mlqtMaxLevel',    0,    ...
             'mlqtLevelSteps',  -1,   ...
			 'CCRho',           @(x) ones(size(x,1),1),...
             'faceConstraints', {{}}, ...
             'FCFactor',        1,  ...
             'interpolateFC',   false, ...
             'circleFactor',    0.6,  ...  
             'protLayer',       false, ...
             'protD',           {{@(p) ones(size(p,1),1)*norm(celldim)/10}},...
             'polyBdr',         zeros(0,2),...
             'useMrstPebi',     false);

opt = merge_options(opt, varargin{:});
circleFactor = opt.circleFactor;

% Set grid sizes
CCGridSize     = min(celldim)*opt.CCFactor;
FCGridSize     = min(celldim)*opt.FCFactor;
mlqtMaxLevel   = opt.mlqtMaxLevel;
mlqtLevelSteps = opt.mlqtLevelSteps;
CCRho          = @(x) CCGridSize*opt.CCRho(x);

% Test input
assert(numel(pdims)==2);
assert(all(pdims>0 ));
assert(CCGridSize>0);
assert(mlqtMaxLevel>=0);
assert(FCGridSize>0);
assert(0.5<circleFactor && circleFactor<1);

if ~isempty(opt.cellConstraints)
    if (numel(opt.interpolateCC) == 1)
        opt.interpolateCC = repmat(opt.interpolateCC, numel(opt.cellConstraints),1);
    end
    assert(numel(opt.interpolateCC)==numel(opt.cellConstraints));
    
    if numel(opt.protD) == 1
        opt.protD = repmat(opt.protD,numel(opt.cellConstraints),1);
    end
    assert(numel(opt.protD) == numel(opt.cellConstraints));
end

if ~isempty(opt.faceConstraints)
    if (numel(opt.interpolateFC) == 1)
        opt.interpolateFC = repmat(opt.interpolateFC, numel(opt.faceConstraints),1);
    end
    assert(numel(opt.interpolateFC)==numel(opt.faceConstraints));
end


if ~all(celldim > 0)
   error('CELLDIM must be positive');
end
if numel(celldim)~=2
  error('CELLDIM must have 2 elements')
end

% Split face constraints and cell constraints
[faceConstraints, fCut, fcCut, IC] = splitAtInt2D(opt.faceConstraints, opt.cellConstraints);
interpFL = opt.interpolateFC(IC);

[cellConstraints,  cCut, cfCut, IC] = splitAtInt2D(opt.cellConstraints, opt.faceConstraints);
interpWP = opt.interpolateCC(IC);
protD    = opt.protD(IC);

% find vertical cell constraints
nw        = cellfun(@numel, opt.cellConstraints)/2;
vW        = nw==1;
cellConstraints = [cellConstraints, opt.cellConstraints(vW)];
cCut      = [cCut; zeros(sum(vW), 1)];
cfCut     = [cfCut; zeros(sum(vW), 1)];
interpWP  = [interpWP; opt.interpolateCC(vW)];
protD     = [protD; opt.protD(vW)];

% Create cell constraint sites
bisectPnt = (FCGridSize.^2 - (circleFactor*FCGridSize).^2 ...
            + (circleFactor*FCGridSize).^2) ./(2*FCGridSize);
faultOffset = sqrt((circleFactor*FCGridSize).^2 - bisectPnt.^2);
sePtn = [cfCut==2|cfCut==3, cfCut==1|cfCut==3];
sePtn = (1.0+faultOffset/CCGridSize)*sePtn;

[CCPts, CCGs,protPts,pGs] = ...
    lineSites2D(cellConstraints, CCGridSize,...
                         'sePtn',         sePtn, ...
                         'cCut',          cCut,...
                         'protLayer',     opt.protLayer,...
                         'protD',         protD, ...
                         'CCRho',         CCRho, ...
                         'interpolateCC', interpWP);

% Create fault points
F = surfaceSites2D(faceConstraints, FCGridSize, ...
                          'circleFactor',  circleFactor, ...
                          'fCut',          fCut, ...
                          'fcCut',         fcCut, ...
                          'interpolateFC', interpFL);

% Create reservoir grid
polyBdr = opt.polyBdr;
[k,n]   = size(polyBdr);
if k==0 % No polygon is given. Assume rectangular box given by pdims
	dx = pdims(1)/ceil(pdims(1)/celldim(1));
	dy = pdims(2)/ceil(pdims(2)/celldim(2));
	vx = 0:dx:pdims(1);
	vy = 0:dy:pdims(2);
    polyBdr = [0, 0; pdims(1), 0; pdims(1), pdims(1); 0, pdims(2)];
elseif k<3
	error('Polygon must have at least 3 edges.');
else
    % A polygon is given, use this and ignore the pdims
	assert(n==2,'polygon boundary is only supported in 2D');
	lDim = [min(polyBdr); max(polyBdr)];
	dx = diff(lDim(:,1))/ceil(diff(lDim(:,1))/celldim(1));
	dy = diff(lDim(:,2))/ceil(diff(lDim(:,2))/celldim(2));
	vx = lDim(1,1):dx:lDim(2,1);
	vy = lDim(1,2):dy:lDim(2,2);
end

[X, Y] = meshgrid(vx, vy);
resPtsInit = [X(:), Y(:)];
if k>=3 % If k < 3 we have a Cartesian box, no need to remove points
	IN         = inpolygon(resPtsInit(:,1),resPtsInit(:,2), polyBdr(:,1), polyBdr(:,2));
    resPtsInit = resPtsInit(IN,:);
end

% Remove tip sites outside domain
if size(F.f.pts, 1) > 0
    innside = inpolygon(F.t.pts(:, 1), F.t.pts(:, 2), polyBdr(:, 1), polyBdr(:, 2));
    F.t.pts = F.t.pts(innside, :);
end

% Refine reservoir grid
if ~isempty(CCPts)
    res = {};
    varArg = {'level', 1, 'maxLev', mlqtMaxLevel, 'distTol', mlqtLevelSteps};
    for i = 1:size(resPtsInit,1)
        res = [res; mlqt(resPtsInit(i,:), CCPts, celldim, varArg{:})];
    end
    resPts = vertcat(res{:, 1});
    %resGridSize = 0.5*[res{:,2}]';
else
    resPts = resPtsInit;
    % resGridSize = repmat(0.5*min(dx,dy),size(resPts,1),1);
end

% Remove conflict points
resPts = removeConflictPoints(resPts, CCPts,  CCGs);
resPts = removeConflictPoints(resPts, protPts,  pGs);
resPts = removeConflictPoints(resPts, F.f.pts,  F.f.Gs);
resPts = removeConflictPoints(resPts, F.c.CC,   F.c.R);

% Create grid
Pts = [F.f.pts; CCPts; protPts; F.t.pts; resPts];

if opt.useMrstPebi
    G = pebi(triangleGrid(Pts));
else
    G = clippedPebi2D(Pts, polyBdr);
end

% Label fault faces.
if ~isempty(F.f.pts)
  N      = G.faces.neighbors + 1; 
  % N == 1 is now a boundary face, so we have to add 1 to the start of 
  % cPos.  We also add empty mapping if there is no constrained faces.
  f2cPos = [1;F.f.cPos; F.f.cPos(end)*ones(size(Pts,1)-size(F.f.pts,1),1)];
  map1   = arrayfun(@colon, f2cPos(N(:,1)),f2cPos(N(:,1)+1)-1,'un',false);
  map2   = arrayfun(@colon, f2cPos(N(:,2)),f2cPos(N(:,2)+1)-1,'un',false);
  G.faces.tag = cellfun(@(c1,c2) numel(intersect(F.f.c(c1),F.f.c(c2)))>1, map1,map2);
else
  G.faces.tag = false(G.faces.num,1);
end

% Label constrained cells
if ~isempty(CCPts)
  G.cells.tag = false(G.cells.num,1);
  
  % Add tag to all cells generated from CCPts
  wellCells = size(F.f.pts,1)+1:size(F.f.pts,1)+size(CCPts,1);
  G.cells.tag(wellCells)= true;
  
  % Add tag to well-fault crossings
  endOfLine = fcCut==1 | fcCut==3;        % Crossing at end of fault
  strOfLine = fcCut==2 | fcCut==3;        % Crossing at start of fault
  fIde      = F.l.fPos([false;endOfLine]);
  fIds      = F.l.fPos(strOfLine);
  fToTag    = [F.l.f(mcolon(fIde - 2,fIde-1)); F.l.f(mcolon(fIds,fIds+1))];

  G.cells.tag(fToTag) = true;
else
  G.cells.tag = false(G.cells.num,1);
end
end
