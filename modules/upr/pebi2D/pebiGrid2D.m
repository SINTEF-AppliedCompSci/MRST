function [G,Pts,F] = pebiGrid2D(resGridSize, pdims, varargin)
% Construct a 2D Pebi grid.
%
% SYNOPSIS:
%   G = pebiGrid2D(resGridSize, pdims)
%   G = pebiGrid2D(...,'Name1',Value1,'Name2',Value2,...)
%
% PARAMETERS
%   resGridSize     - Size of the reservoir grid cells, in units of
%                     meters.
%
%   pdims           - Vector, length 2, [xmax, ymax], of physical size in
%                     units of meters of the computational domain.
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
%   CCFactor  - OPTIONAL.
%                     Default value is 1. This gives the relative grid
%                     size of the cell constraint grid cells compared to
%                     reservoir grid cells. If CCFactor=0.5 the constrained
%                     cells will be about half the size of the reservoir cells.
%
%   CCRefinement  - OPTIONAL
%                     Default value FALSE. Set to true to turn on
%                     refinement around the cell constraints.
%
%   CCEps         - OPTIONAL
%                     Default value 0.25/max(pdims). CCEps set the
%                     refinement transition around cell constraints. The density
%                     function for the reservoir grid is set by
%                     rho~exp(-distance to well / CCEps).
%
%   CCRho         - OPTIONAL
%                     Default value @(x) ones(size(x,1),1). Function gives
%                     the relative distance between cell constraint cells. If
%                     CCRho=0.5 in an area the distance between
%                     constrained-cells will be 0.5*CCFactor*resGridSize
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
%    protD          - OPTIONAL.
%                     Default value CCGridSize/10. Cell array of Functions.
%                     The array should have either one function or one
%                     function for  each well path.
%                     The functions give the distance from cell constrainted  sites
%                     to the protection sites. The function is evaluated along
%                     the well path such that protD(0) is the start of the
%                     well while protD(1) is the end of the well.
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
%   FCRho           - OPTIONAL.
%                     Default value 1. Function that gives the relative
%                     distance between surfaces sites along a path. If
%                     FCRho=0.5 in a area, the surface sites will here
%                     be about 50% closer than in other areas.
%
%   FCRefinement - OPTIONAL
%                     Default value FALSE. Set to true to turn on
%                     refinement around face constraints.
%
%   FCEps        - OPTIONAL
%                     Default value 0.25/max(pdims). FCEps set the
%                     refinement transition around face constraints.
%                     The density function for the reservoir grid is set by
%                     rho~exp(-distance to constraint / FCEps).
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
%   sufFCCond    - OPTIONAL
%                     Default value true. If sufFCCond = false we do not
%                     enforce the sufficient and necessary face constraint condition.
%                     Instead we enforce a less strict condition and remove
%                     any reservoir sites that are closer to the constraint
%                     sites than the constraint grid size. Note that Conformity
%                     is then not guaranteed. You might still set this to
%                     false if you have problems with bad cells at the end
%                     of your faults because the sufficient condition
%                     removes some reservoir points.
%
%   linearize       - OPTIONAL
%                     Default is false. If true, we evaluate the scaled
%                     edge-length function in distmesh2d by interpolating
%                     between values computed at the vertices. There are
%                     usually more edges than vertices and if the
%                     edge-length function is expensive to compute, this
%                     will save significant computational time.
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
%   G  = pebiGrid2D(1/10,[1,1],'cellConstraints',wl,'faceConstraints',fl)
%   cla, plotGrid(G)
%
% SEE ALSO
%   compositePebiGrid2D, pebi, surfaceSites2D.

%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2015-2020 Runar Lie Berge. See COPYRIGHT.TXT for details.
%
% distMesh is used to create the background grid
% (http://persson.berkeley.edu/distmesh/)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}

% Set options
opt = struct('cellConstraints', {{}}, ...
             'CCFactor',        1, ...
             'interpolateCC',   false, ...
             'CCRefinement',    false, ...
			 'CCRho',           @(x) ones(size(x,1),1),...
             'CCEps',           -1,...
             'protLayer',       false, ...
             'protD',           {{@(p) ones(size(p,1),1)*resGridSize/10}},...
             'faceConstraints', {{}}, ...
             'FCFactor',        1, ...
             'interpolateFC',   false, ...
             'circleFactor',    0.6,...
             'FCRefinement',    false, ...
             'FCRho',           @(x) ones(size(x,1),1),...
             'FCEps',           -1,...
			 'polyBdr',         zeros(0,2), ...
			 'sufFCCond',       true,...,
             'linearize',       false,...,
             'useMrstPebi',     false);

opt          = merge_options(opt, varargin{:});
circleFactor = opt.circleFactor;
CCRef        = opt.CCRefinement;
FCRef        = opt.FCRefinement;
CCEps        = opt.CCEps;
FCEps        = opt.FCEps;
% Set grid sizes
CCFactor   = opt.CCFactor;
FCFactor   = opt.FCFactor;
CCGridSize = resGridSize*CCFactor;
CCRho      = @(x) CCGridSize*opt.CCRho(x);
FCGridSize = resGridSize*FCFactor;
FCRho      = opt.FCRho;

% If polygon boundary is give, use this
polyBdr = opt.polyBdr;
[sizePolyBdr,l] = size(polyBdr);
if sizePolyBdr >= 3
    pdims = max(polyBdr) - min(polyBdr);
elseif sizePolyBdr==2 || sizePolyBdr == 1
	error('Polygon must have at least 3 edges.');
end

if CCEps<0
    CCEps = 0.25 * max(pdims);
end
if FCEps<0
    FCEps = 0.25 * max(pdims);
end

% Test input
assert(resGridSize>0);
assert(numel(pdims)==2);
assert(all(pdims>0 ));
assert(CCGridSize>0);
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

% Split face constraints and cell constraints
[faceConstraints, fCut, fcCut, IC] = splitAtInt2D(opt.faceConstraints, opt.cellConstraints);
interpFL = opt.interpolateFC(IC);

[cellConstraints,  cCut, cfCut, IC] = splitAtInt2D(opt.cellConstraints, opt.faceConstraints);
interpWP = opt.interpolateCC(IC);
protD    = opt.protD(IC);

% find vertical cell constraints
nw        = cellfun(@numel, opt.cellConstraints)/2;
vW        = nw==1;
cellConstraints = [cellConstraints,opt.cellConstraints(vW)];
cCut      = [cCut;zeros(sum(vW),1)];
cfCut     = [cfCut; zeros(sum(vW),1)];
interpWP  = [interpWP; opt.interpolateCC(vW)];
protD     = [protD; opt.protD(vW)];

% Create cell constraint sites
sePtn = [cfCut==2|cfCut==3, cfCut==1|cfCut==3];
if CCRef
  fLen = CCGridSize*FCFactor*1.2;  
else
  fLen = FCGridSize;
end
bisectPnt = (fLen.^2 - (circleFactor*fLen).^2 ...
            + (circleFactor*fLen).^2) ./(2*fLen);
faultOffset = sqrt((circleFactor*fLen).^2 - bisectPnt.^2);
sePtn = (1.0+faultOffset/CCGridSize)*sePtn;

[wellPts, wGs, protPts, pGs] = ...
    lineSites2D(cellConstraints, CCGridSize, ...
                         'sePtn',         sePtn,...
                         'cCut',          cCut, ...
                         'protLayer',     opt.protLayer,...
                         'protD',         protD, ...
                         'CCRho',         CCRho, ...
                         'interpolateCC', interpWP);


% create distance functions
if CCRef && ~isempty(wellPts)
    hresw  = @(x) min((ones(size(x,1),1)/CCFactor), ...
        1.2*exp(minPdist2(x,wellPts)/CCEps));
    hfault = @(x) CCGridSize*FCFactor*hresw(x).*FCRho(x);
else
    hresw  = @(p) constFunc(p)/CCFactor;
    hfault = @(p) FCGridSize*FCRho(p);
end

% Create surface sites
F = surfaceSites2D(faceConstraints, FCGridSize,...
                          'circleFactor',  circleFactor,...
                          'fCut',          fCut, ...
                          'fcCut',         fcCut, ...
                          'interpolateFC', interpFL, ...
                          'distFun',       hfault);

if FCRef && ~isempty(F.f.pts)
  hresf = @(x) min((ones(size(x,1),1)/FCFactor), ...
                    1.2*exp(minPdist2(x,F.f.pts)/FCEps));
else
  hresf = @(p) constFunc(p)/CCFactor;
end

% Create reservoir grid points
% set domain function
if sizePolyBdr==0
	x = pdims;
	rectangle = [0,0; x(1),x(2)];
	fd = @(p,varargin) drectangle(p, 0, x(1), 0, x(2));
	corners = [0,0; 0,x(2); x(1),0; x(1),x(2)];
	vararg  = [];
    polyBdr = [0, 0; pdims(1), 0; pdims(1), pdims(2); 0, pdims(2)];
else
    assert(l==2,'polygon boundary is only supported in 2D');
	rectangle = [min(polyBdr); max(polyBdr)];
	corners   = polyBdr;
	fd        = @dpoly;
	vararg    = [polyBdr; polyBdr(1,:)];
end

% Remove tip sites outside domain
if size(F.t.pts, 1) > 0
    innside = inpolygon(F.t.pts(:, 1), F.t.pts(:, 2), polyBdr(:, 1), polyBdr(:, 2));
    F.t.pts = F.t.pts(innside, :);
end    

if FCRef && CCRef
    ds   = min(CCGridSize,FCGridSize);
    hres = @(x,varargin) min(hresf(p), hresw(p));
elseif FCRef
    ds   = FCGridSize;
    hres = @(p,varargin) hresf(p);
else 
    ds   = CCGridSize;
    hres = @(p, varargin) hresw(p);
end
fixedPts = [F.f.pts; wellPts; protPts; F.t.pts; corners];  
[Pts,~,sorting] = distmesh2d(fd, hres, ds, rectangle, fixedPts, opt.linearize, vararg);

% Distmesh change the order of all sites. We undo this sorting.
isFault = false(max(sorting),1); isFault(1:size(F.f.pts,1)) = true;
isFault = isFault(sorting);
[~,If]  = sort(sorting(isFault));
isWell  = false(max(sorting),1); isWell(size(F.f.pts,1)+(1:size(wellPts,1)))= true;
isWell  = isWell(sorting);
[~,Iw]  = sort(sorting(isWell));
isRes   = ~isFault & ~isWell;

fPts = Pts(isFault,:);
fPts = fPts(If,:);
wPts = Pts(isWell,:);
wPts = wPts(Iw,:);

if opt.sufFCCond
	Pts = surfaceSufCond2D(Pts(isRes,:),F);
else
	Pts = removeConflictPoints(Pts(isRes,:),F.f.pts, F.f.Gs);
end
Pts  = [fPts; wPts; Pts];
% Create grid
if opt.useMrstPebi
	t    = delaunay(Pts);
	% Fix boundary
	pmid = (Pts(t(:,1),:)+Pts(t(:,2),:)+Pts(t(:,3),:))/3;% Compute centroids
	t    = t(fd(pmid,vararg)<-0.001*CCFactor,:);   % Keep interior triangles

	G = triangleGrid(Pts, t);
	G = pebi(G);
else
    G = clippedPebi2D(Pts, polyBdr);
end
% label face constraint faces.
if ~isempty(F.f.pts)
  N      = G.faces.neighbors + 1; 
  % N == 1 is now a boundary constraint, so we have to add 1 to the start of 
  % cPos.  We also add empty mapping for no surface pts.
  f2cPos = [1;F.f.cPos; F.f.cPos(end)*ones(size(Pts,1)-size(F.f.pts,1),1)];
  map1   = arrayfun(@colon, f2cPos(N(:,1)),f2cPos(N(:,1)+1)-1,'un',false);
  map2   = arrayfun(@colon, f2cPos(N(:,2)),f2cPos(N(:,2)+1)-1,'un',false);
  G.faces.tag = cellfun(@(c1,c2) numel(intersect(F.f.c(c1),F.f.c(c2)))>1, map1,map2);
else
  G.faces.tag = false(G.faces.num,1);
end

% Label cell constraint cells
if ~isempty(wellPts)
  G.cells.tag = false(G.cells.num,1);
  % Add tag to all cells generated from wellPts
  wellCells = size(F.f.pts,1)+1:size(F.f.pts,1)+size(wellPts,1);
  G.cells.tag(wellCells)= true;
  
  % Add tag to cell and face constraint crossings
  endOfLine = fcCut==1 | fcCut==3;        % Crossing at end of face constraint
  strOfLine = fcCut==2 | fcCut==3;        % Crossing at start of face constraint
  fIde      = F.l.fPos([false;endOfLine]);
  fIds      = F.l.fPos(strOfLine);
  fToTag    = [F.l.f(mcolon(fIde - 2,fIde-1)); F.l.f(mcolon(fIds,fIds+1))];

  G.cells.tag(fToTag) = true;
else
  G.cells.tag = false(G.cells.num,1);
end
end
