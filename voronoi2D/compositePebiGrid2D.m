function [G,Pts,F] = compositePebiGrid2D(celldim, pdims, varargin)
% Construct a 2D composite Pebi grid. A cartesian background grid that is 
% refined around faults and wells.
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
%   wellLines         - OPTIONAL.
%                       Default value empty. A struct of vectors. Each 
%                       vector, size nw x 2, is the coordinates of a 
%                       well-trace. The well is assumed to be linear 
%                       between the coorinates. If the vector only contains 
%                       one coordinate, the well is treated as a point well.
%
%   wellGridFactor    - OPTIONAL.
%                       Default value is 0.5. This gives the relative grid
%                       size of the well grid cells compared to reservoir 
%                       grid cells. If wellGridFactor=0.5 the well cells 
%                       will be about half the size of the reservoir cells.
%
%   interpolWP        - OPTIONAL.
%                       Default value is a boolean false value, but the
%                       user can supply an individual value per well path.
%                       If false, each segment in the corresponding well
%                       curve will be represented by at least one cell. If
%                       true, the routine will interpolate along the
%                       curve, which means that cell centers will not
%                       necessarily fall exactly on the prescribed curve.
%
%   mlqtMaxLevel      - OPTIONAL.
%                       Default value 0. Number of refinement steps around 
%                       wells. 
%
%   mlqtLevelSteps    - OPTIONAL.  
%                       Default value -1. Vector of length mlqtMaxLevel 
%                       which specifies the radius of each refinement level.
%                       The default value -1 calls the default level step
%                       in the mlqt function.
%
%   wellRho           - OPTIONAL
%                       Default value @(x) ones(size(x,1),1). Function gives
%                       the relative distance between well points. If
%                       wellRho=0.5 in an area, the distance between
%                       well-cells will be
%                       0.5*wellGridFactor*min(resGridSize)
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
%    protD            - OPTIONAL.
%                       Default value wellGridSize/10. Cell array of Functions.
%                       The array should have either one function or one 
%                       function for  each well path.
%                       The functions give the distance from well sites 
%                       to the protection sites. The function is evaluated along 
%                       the well path such that protD(0) is the start of the 
%                       well while protD(1) is the end of the well.
%
%   faultLines        - OPTIONAL
%                       Default value empty. A struct of vectors.  Each 
%                       vector, size nf x 2, is the coordinates of a 
%                       fault-trace. The fault is assumed to be linear 
%                       between the coorinates
%
%   faultGridFactor   - OPTIONAL.
%                       Default value is 0.5. This gives the relative grid
%                       size of the fault grid cells compared to reservoir 
%                       grid cells. If faultGridFactor=0.5, the fault cells 
%                       will be about half the size of the reservoir cells.
%
%   interpolFL        - OPTIONAL.
%                       Default value is a boolean false value, but the
%                       user can supply an individual value per well path.
%                       If false, each segment in the corresponding fault
%                       curve will be represented by at least one cell. If
%                       true, the routine will interpolate along the
%                       curve, which means that cell edges will not
%                       necessarily fall exactly on the prescribed curve.
%
%   circleFactor      - OPTIONAL.
%                       Default value 0.6.  Valid values are between 0.5 
%                       and 1. circleFactor controll the size of the 
%                       circles used to create the fault grid points. The 
%                       circleFactor is the ratio between the radius, 
%                       and distace between the circles. A small value will
%                       place the fault points close the the faults, while
%                       a large value will place the far from the faults.
%
%   polyBdr           - OPTIONAL 
%                       Default value []. plyBdr is a array of size [k,2].
%                       if k>=3 polyBdr gives the vertices of the reservoir
%                       boundary. For this domain:
%                                  .(x1,y1)
%                                 / \ 
%                         (x3,y3).---.(x2,y2)
%                       polyBdr would be [x1,y1;x2,y2;x3,y3]. The set of
%                       vertices must run clockwise or counterclockwise.
%                       Note that if k>=3 the innput pdims will have no
%                       effect.
%
% RETURNS:
%   G                - Valid grid definition.  
%                        Contains extra fields:
%                          - G.cells.tag is TRUE for all well cells.
%                          - G.faces.tag is TRUE for all fault edges.
%   Pts              - Array [G.cells.num x 3] of the Voronoi sites.
%   F                - Struct with elements as returned from 
%                      surfaceSites2D
%
% EXAMPLE:
%   fl = {[0.2,0.2;0.8,0.8]};
%   wl = {[0.2,0.8;0.8,0.2]};
%   G  = compositePebiGrid2D([1/10,1/10],[1,1],'wellLines',wl,'faultLines',fl)
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
opt = struct('wellLines',       {{}}, ...
             'wellGridFactor',  1,  ...
             'interpolWP',      false, ...
             'mlqtMaxLevel',    0,    ...
             'mlqtLevelSteps',  -1,   ...
			 'wellRho',         @(x) ones(size(x,1),1),...
             'faultLines',      {{}}, ...
             'faultGridFactor', 1,  ...
             'interpolFL',      false, ...
             'circleFactor',    0.6,  ...  
             'protLayer',       false, ...
             'protD',           {{@(p) ones(size(p,1),1)*norm(celldim)/10}},...
             'polyBdr',         zeros(0,2),...
             'useMrstPebi',     false);

opt = merge_options(opt, varargin{:});
circleFactor = opt.circleFactor;

% Set grid sizes
wellGridSize   = min(celldim)*opt.wellGridFactor;
faultGridSize  = min(celldim)*opt.faultGridFactor;
mlqtMaxLevel   = opt.mlqtMaxLevel;
mlqtLevelSteps = opt.mlqtLevelSteps;
wellRho        = @(x) wellGridSize*opt.wellRho(x);

% Test input
assert(numel(pdims)==2);
assert(all(pdims>0 ));
assert(wellGridSize>0);
assert(mlqtMaxLevel>=0);
assert(faultGridSize>0);
assert(0.5<circleFactor && circleFactor<1);

if ~isempty(opt.wellLines)
    if (numel(opt.interpolWP) == 1)
        opt.interpolWP = repmat(opt.interpolWP, numel(opt.wellLines),1);
    end
    assert(numel(opt.interpolWP)==numel(opt.wellLines));
    
    if numel(opt.protD) == 1
        opt.protD = repmat(opt.protD,numel(opt.wellLines),1);
    end
    assert(numel(opt.protD) == numel(opt.wellLines));
end

if ~isempty(opt.faultLines)
    if (numel(opt.interpolFL) == 1)
        opt.interpolFL = repmat(opt.interpolFL, numel(opt.faultLines),1);
    end
    assert(numel(opt.interpolFL)==numel(opt.faultLines));
end


if ~all(celldim > 0)
   error('CELLDIM must be positive');
end
if numel(celldim)~=2
  error('CELLDIM must have 2 elements')
end

% Split faults and wells paths
[faultLines, fCut, fwCut, IC] = splitAtInt(opt.faultLines, opt.wellLines);
interpFL = opt.interpolFL(IC);

[wellLines,  wCut, wfCut, IC] = splitAtInt(opt.wellLines, opt.faultLines);
interpWP = opt.interpolWP(IC);
protD    = opt.protD(IC);

% find vertical wells
nw        = cellfun(@numel, opt.wellLines)/2;
vW        = nw==1;
wellLines = [wellLines,opt.wellLines(vW)];
wCut      = [wCut;zeros(sum(vW),1)];
wfCut     = [wfCut; zeros(sum(vW),1)];
interpWP  = [interpWP; opt.interpolWP(vW)];
protD     = [protD; opt.protD(vW)];

% Create well points
bisectPnt = (faultGridSize.^2 - (circleFactor*faultGridSize).^2 ...
            + (circleFactor*faultGridSize).^2) ./(2*faultGridSize);
faultOffset = sqrt((circleFactor*faultGridSize).^2 - bisectPnt.^2);
sePtn = [wfCut==2|wfCut==3, wfCut==1|wfCut==3];
sePtn = (1.0+faultOffset/wellGridSize)*sePtn;

[wellPts, wGs,protPts,pGs] = ...
    lineSites2D(wellLines, wellGridSize,...
                         'sePtn',        sePtn, ...
                         'wCut',         wCut,...
                         'protLayer',    opt.protLayer,...
                         'protD',        protD, ...
                         'wellRho',      wellRho, ...
                         'interpolWP',   interpWP);

% Create fault points
F = surfaceSites2D(faultLines, faultGridSize, ...
                          'circleFactor', circleFactor, ...
                          'fCut',         fCut, ...
                          'fwCut',        fwCut, ...
                          'interpolFL',   interpFL);

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
if ~isempty(wellPts)
    res = {};
    varArg = {'level', 1, 'maxLev', mlqtMaxLevel, 'distTol', mlqtLevelSteps};
    for i = 1:size(resPtsInit,1)
        res = [res; mlqt(resPtsInit(i,:), wellPts, celldim, varArg{:})];
    end
    resPts = vertcat(res{:, 1});
    %resGridSize = 0.5*[res{:,2}]';
else
    resPts = resPtsInit;
    % resGridSize = repmat(0.5*min(dx,dy),size(resPts,1),1);
end

% Remove Conflic Points
resPts = removeConflictPoints2(resPts, wellPts,  wGs);
resPts = removeConflictPoints2(resPts, protPts,  pGs);
resPts = removeConflictPoints2(resPts, F.f.pts,  F.f.Gs);
resPts = removeConflictPoints2(resPts, F.c.CC,   F.c.R);

% Create Grid
Pts = [F.f.pts; wellPts; protPts; F.t.pts; resPts];

if opt.useMrstPebi
    G = pebi(triangleGrid(Pts));
else
    G = clippedPebi2D(Pts, polyBdr);
end

% label fault faces.
if ~isempty(F.f.pts)
  N      = G.faces.neighbors + 1; 
  % N == 1 is now a boundary fault, so we have to add 1 to the start of 
  % cPos.  We also add empty mapping for no fault pts.
  f2cPos = [1;F.f.cPos; F.f.cPos(end)*ones(size(Pts,1)-size(F.f.pts,1),1)];
  map1   = arrayfun(@colon, f2cPos(N(:,1)),f2cPos(N(:,1)+1)-1,'un',false);
  map2   = arrayfun(@colon, f2cPos(N(:,2)),f2cPos(N(:,2)+1)-1,'un',false);
  G.faces.tag = cellfun(@(c1,c2) numel(intersect(F.f.c(c1),F.f.c(c2)))>1, map1,map2);
else
  G.faces.tag = false(G.faces.num,1);
end

%Label well cells
if ~isempty(wellPts)
  G.cells.tag = false(G.cells.num,1);
  
  % Add tag to all cells generated from wellPts
  wellCells = size(F.f.pts,1)+1:size(F.f.pts,1)+size(wellPts,1);
  G.cells.tag(wellCells)= true;
  
  % Add tag to well-fault crossings
  endOfLine = fwCut==1 | fwCut==3;        % Crossing at end of fault
  strOfLine = fwCut==2 | fwCut==3;        % Crossing at start of fault
  fIde      = F.l.fPos([false;endOfLine]);
  fIds      = F.l.fPos(strOfLine);
  fToTag    = [F.l.f(mcolon(fIde - 2,fIde-1)); F.l.f(mcolon(fIds,fIds+1))];

  G.cells.tag(fToTag) = true;
else
  G.cells.tag = false(G.cells.num,1);
end
end
