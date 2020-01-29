function [G,Pts,F] = pebiGrid(resGridSize, pdims, varargin)
% Construct a 2D Pebi grid.
%
% SYNOPSIS:
%   G = pebiGrid(resGridSize, pdims)
%   G = pebiGrid(...,'Name1',Value1,'Name2',Value2,...)
%
% PARAMETERS
%   resGridSize     - Size of the reservoir grid cells, in units of
%                     meters. 
%
%   pdims           - Vector, length 2, [xmax, ymax], of physical size in
%                     units of meters of the computational domain. 
%
%   wellLines       - OPTIONAL.
%                     Default value empty. A struct of vectors. Each 
%                     vector, size nw x 2, is the coordinates of a 
%                     well-trace. The well is assumed to be linear 
%                     between the coorinates. If the vector only contains 
%                     one coordinate, the well is treated as a point well.
%
%   interpolWP     - OPTIONAL.
%                     Default value is a boolean false value, but the user
%                     can supply an individual value per well path. If
%                     false, each segment in the corresponding well curve
%                     will be represented by at least one cell. If true,
%                     the routine will interpolate along the curve, which
%                     means that cell centers will not necessarily fall
%                     exactly on the prescribed curve.
%
%   wellGridFactor  - OPTIONAL.
%                     Default value is 1. This gives the relative grid
%                     size of the well grid cells compared to reservoir 
%                     grid cells. If wellGridFactor=0.5 the well cells 
%                     will be about half the size of the reservoir cells.
%
%   wellRefinement  - OPTIONAL
%                     Default value FALSE. Set to true to turn on
%                     refinement around wells.
%
%   wellEps         - OPTIONAL
%                     Default value 0.25/max(pdims). wellEps set the
%                     refinement transition around wells. The density
%                     function for the reservoir grid is set by
%                     rho~exp(-distance to well / wellEps).
%
%   wellRho         - OPTIONAL
%                     Default value @(x) ones(size(x,1),1). Function gives
%                     the relative distance between well points. If
%                     wellRho=0.5 in an area the distance between
%                     well-cells will be 0.5*wellGridFactor*resGridSize
%
%   protLayer       - OPTIONAL.
%                     Default set to false. If set to true a protection layer
%                     is added on both sides of the well
%                     .          .             .  Protection Layer
%                                                 protD(dist between points)
%                     .----------.-------------.  Well path
%
%                     .          .             .  Protection Layer
%     
%    protD          - OPTIONAL.
%                     Default value wellGridSize/10. Cell array of Functions.
%                     The array should have either one function or one 
%                     function for  each well path.
%                     The functions give the distance from well sites 
%                     to the protection sites. The function is evaluated along 
%                     the well path such that protD(0) is the start of the 
%                     well while protD(1) is the end of the well.
%
%   faultLines      - OPTIONAL
%                     Default value empty. A struct of vectors.  Each 
%                     vector, size nf x 2, is the coordinates of a 
%                     fault-trace. The fault is assumed to be linear 
%                     between the coorinates
%
%   faultGridFactor - OPTIONAL.
%                     Default value is 0.5. This gives the relative grid
%                     size of the fault grid cells compared to reservoir 
%                     grid cells. If faultGridFactor=0.5 the fault cells 
%                     will be about half the size of the reservoir cells.
%
%   interpolFL      - OPTIONAL.
%                     Default value is a boolean false value, but the
%                     user can supply an individual value per well path.
%                     If false, each segment in the corresponding fault
%                     curve will be represented by at least one cell. If
%                     true, the routine will interpolate along the
%                     curve, which means that cell edges will not
%                     necessarily fall exactly on the prescribed curve.
%
%   circleFactor    - OPTIONAL.
%                     Default value 0.6.  Valid values are between 0.5 
%                     and 1. circleFactor controll the size of the 
%                     circles used to create the fault grid points. The 
%                     circleFactor is the ratio between the radius 
%                     and distance between the circles. A small value will
%                     place the fault points close the the faults, while
%                     a large value will place the far from the faults.
%
%   faultRho        - OPTIONAL.
%                     Default value 1. Function that gives the relative
%                     distance between fault sites along a path. If
%                     faultRho=0.5 in a area, the fault sites will here
%                     be about 50% closer than in other areas.
%
%   faultRefinement - OPTIONAL
%                     Default value FALSE. Set to true to turn on
%                     refinement around faults.
%
%   faultEps        - OPTIONAL
%                     Default value 0.25/max(pdims). faultEps set the
%                     refinement transition around faults. The density
%                     function for the reservoir grid is set by
%                     rho~exp(-distance to fault / faultEps).
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
%   sufFaultCond    - OPTIONAL 
%                     Default value true. If sufFaultCond = false we do not
%                     enforce the sufficient and necessary fault condition.
%                     Instead we enforce a less strict condition and remove
%                     any reservoir sites that are closer to the fault
%                     sites than the fault grid size. Note that Conformity
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
%                        - G.cells.tag is TRUE for all well cells.
%                        - G.faces.tag is TRUE for all fault edges.
%   Pts            - Array [G.cells.num x 3] of the Voronoi sites.
%   F              - Struct with elements as returned from 
%                    surfaceSites2D
%
% EXAMPLE:
%   fl = {[0.2,0.2;0.8,0.8]};
%   wl = {[0.2,0.8;0.8,0.2]};
%   G  = pebiGrid(1/10,[1,1],'wellLines',wl,'faultLines',fl)
%   cla, plotGrid(G)
%
% SEE ALSO
%   compositePebiGrid2D, pebi, surfaceSites2D, lineSites2D.

%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2015-2020 Runar Lie Berge. See COPYRIGHT.TXT for details.
%
% distMesh is used to create the background grid 
% (http://persson.berkeley.edu/distmesh/)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}  

% Set options
opt = struct('wellLines',       {{}}, ...
             'wellGridFactor',  1, ...
             'interpolWP',      false, ...
             'wellRefinement',  false, ...
			 'wellRho',         @(x) ones(size(x,1),1),...
             'wellEps',         -1,...
             'protLayer',       false, ...
             'protD',           {{@(p) ones(size(p,1),1)*resGridSize/10}},...
             'faultLines',      {{}}, ...
             'faultGridFactor', 1, ...
             'interpolFL',      false, ...
             'circleFactor',    0.6,...
             'faultRefinement', false, ...
             'faultRho',        @(x) ones(size(x,1),1),...
             'faultEps',        -1,...
			 'polyBdr',         zeros(0,2), ...
			 'sufFaultCond',    true,...,
             'linearize',       false,...,
             'useMrstPebi',     false);

opt             = merge_options(opt, varargin{:});
circleFactor    = opt.circleFactor;
wellRef         = opt.wellRefinement;
faultRef        = opt.faultRefinement;
wellEps         = opt.wellEps;
faultEps        = opt.faultEps;
% Set grid sizes
wellGridFactor  = opt.wellGridFactor;
faultGridFactor = opt.faultGridFactor;
wellGridSize    = resGridSize*wellGridFactor;
wellRho         = @(x) wellGridSize*opt.wellRho(x);
faultGridSize   = resGridSize*faultGridFactor;
faultRho = opt.faultRho;

% If polygon boundary is give, use this
polyBdr = opt.polyBdr;
[sizePolyBdr,l] = size(polyBdr);
if sizePolyBdr >= 3
    pdims = max(polyBdr) - min(polyBdr);
elseif sizePolyBdr==2 || sizePolyBdr == 1
	error('Polygon must have at least 3 edges.');
end

if wellEps<0
    wellEps = 0.25 * max(pdims);
end
if faultEps<0
    faultEps = 0.25 * max(pdims);
end

% Test input
assert(resGridSize>0);
assert(numel(pdims)==2);
assert(all(pdims>0 ));
assert(wellGridSize>0);
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

% Split faults and wells paths
[faultLines, fCut, fwCut, IC] = splitAtInt2D(opt.faultLines, opt.wellLines);
interpFL = opt.interpolFL(IC);

[wellLines,  wCut, wfCut, IC] = splitAtInt2D(opt.wellLines, opt.faultLines);
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

% Create well Points
sePtn = [wfCut==2|wfCut==3, wfCut==1|wfCut==3];
if wellRef
  fLen = wellGridSize*faultGridFactor*1.2;  
else
  fLen = faultGridSize;
end
bisectPnt = (fLen.^2 - (circleFactor*fLen).^2 ...
            + (circleFactor*fLen).^2) ./(2*fLen);
faultOffset = sqrt((circleFactor*fLen).^2 - bisectPnt.^2);
sePtn = (1.0+faultOffset/wellGridSize)*sePtn;

[wellPts, wGs, protPts, pGs] = ...
    lineSites2D(wellLines, wellGridSize, ...
                         'sePtn',        sePtn,...
                         'wCut',         wCut, ...
                         'protLayer',    opt.protLayer,...
                         'protD',        protD, ...
                         'wellRho',      wellRho, ...
                         'interpolWP',   interpWP);


% create distance functions
if wellRef && ~isempty(wellPts)
    hresw  = @(x) min((ones(size(x,1),1)/wellGridFactor), ...
        1.2*exp(minPdist2(x,wellPts)/wellEps));
    hfault = @(x) wellGridSize*faultGridFactor*hresw(x).*faultRho(x);
else
    hresw  = @(p) constFunc(p)/wellGridFactor;
    hfault = @(p) faultGridSize*faultRho(p);
end

% Create fault points
F = surfaceSites2D(faultLines, faultGridSize,...
                          'circleFactor',  circleFactor,...
                          'fCut',          fCut, ...
                          'fwCut',         fwCut, ...
                          'interpolFL',   interpFL, ...
                          'distFun', hfault);

if faultRef && ~isempty(F.f.pts)
  hresf = @(x) min((ones(size(x,1),1)/faultGridFactor), ...
                    1.2*exp(minPdist2(x,F.f.pts)/faultEps));
else
  hresf = @(p) constFunc(p)/wellGridFactor;
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

if faultRef && wellRef
    ds   = min(wellGridSize,faultGridSize);
    hres = @(x,varargin) min(hresf(p), hresw(p));
elseif faultRef
    ds   = faultGridSize;
    hres = @(p,varargin) hresf(p);
else 
    ds   = wellGridSize;
    hres = @(p, varargin) hresw(p);
end
fixedPts = [F.f.pts; wellPts; protPts; F.t.pts; corners];  
[Pts,~,sorting] = distmesh2d(fd, hres, ds, rectangle, fixedPts, opt.linearize, vararg);

% Distmesh change the order of all points. We undo this sorting.
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

if opt.sufFaultCond
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
	t    = t(fd(pmid,vararg)<-0.001*wellGridFactor,:);   % Keep interior triangles

	G = triangleGrid(Pts, t);
	G = pebi(G);
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
