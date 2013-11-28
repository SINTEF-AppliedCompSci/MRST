function [Gt, rock, rock2D, bcIxVE] = makeSleipnerVEmodelSimple(usemex)
%Make an VE model based upon the Sleipner data set from ieaghg.org
%
% SYNOPSIS:
%  [G, Gt, bcIx, bcIxVE, rock, rock2D] = makeSleipnerVEmodel()
%
% PARAMETERS:
%   G      - Data structure for 3D grid
%   Gt     - Data structure for topsurface grid
%   rock   - Data structure for 3D rock parameters
%   rock2D - Data structure for rock parameters for topsurface grid
%   bcIxVE - Index for pressure boundary conditions in topsurface grid
%   usemex - Flag: if true, use C-accelerated routines for processing
%            Eclipse input and computing geometry
%
% DESCRIPTION:
%
% SEE ALSO:
%   runSleipner

if nargin<1, usemex = true; end

try
   disp(' -> Reading SleipnerSimple.mat');
   datadir = fullfile(VEROOTDIR, 'data', 'mat');
   load(fullfile(datadir,'SleipnerSimple'));
   return;
catch %#ok<*CTCH>
   disp('    Data set has not yet been constructed.');
end

%% Read data
try
   sdir = fullfile('data', 'sleipner');
   disp([' -> Reading data from: ' sdir]);
   grdecl = readGRDECL(fullfile(VEROOTDIR, sdir, 'SLEIPNER.DATA'));
catch
   fprintf(1, '    Reading failed, please dowload data manually following');
   fprintf(1, ' instructions\n    in "%s"\n', fullfile(sdir,'README'));
   return;
end

%% Process 3D grid and compute geometry
% First, we map from left-hand to right-hand coordinate system. 
disp(' -> Processing grid');
lines = reshape(grdecl.COORD,6,[]);
lines([2 5],:) = -lines([2 5],:);
grdecl.COORD = lines(:); clear lines

% Then, we remove the bottom and top layers that contain shale
%grdecl.ACTNUM(grdecl.PERMX<200) = 0;
if(usemex)
  G=mprocessGRDECL(grdecl);
else
  G=processGRDECL(grdecl); 
end
[ijk{1:3}]=ind2sub(G.cartDims,G.cells.indexMap)
kmin=double(min(ijk{3}(grdecl.PERMX>200)));
kmax=double(max(ijk{3}(grdecl.PERMX>200)));
grdecl=cutGrdecl(grdecl,[1,grdecl.cartDims(1);1,grdecl.cartDims(2);kmin,kmax]);
permx=mean(grdecl.PERMX);poro=mean(grdecl.PORO);
grdecl=coarseGrdecl(grdecl,[1 1,grdecl.cartDims(3)]);
grdecl.PERMX=permx*ones(size(grdecl.ACTNUM));
grdecl.PERMY=permx*ones(size(grdecl.ACTNUM));
grdecl.PORO=poro*ones(size(grdecl.ACTNUM));
%grdecl = coarseGrdecl(grdecl,[1 1,
% Next, we process the grid and compute geometry, possibly using
% C-accelerated routines
if usemex,
   mlist = mrstModule;
   mrstModule add 'mex/opm_gridprocessing'
   mrstModule add 'mex/libgeometry'
   G = mprocessGRDECL(grdecl);
   G = mcomputeGeometry(G);
   mrstModule('reset', mlist{:});
else
   G = processGRDECL(grdecl);
   G = computeGeometry(G);
end

% Adding tags needed by topSurfaceGrid
G.cells.faces = [G.cells.faces, repmat((1:6).', [G.cells.num, 1])];


%% Construct petrophysical model
rock = grdecl2Rock(grdecl, G.cells.indexMap);
rock.perm = convertFrom(rock.perm, milli*darcy);
clear grdecl


%% Construct top-surface grid
disp(' -> Constructing top-surface grid');
[Gt, G] = topSurfaceGrid(G);
rock2D  = averageRock(rock, Gt);


%% Find pressure boundary
% Setting boundary conditions is unfortunately a manual process and may
% require some fiddling with indices, as shown in the code below. Here, we
% need to find all outer vertical faces
i = any(Gt.faces.neighbors==0, 2);  % find all outer faces
I = i(Gt.cells.faces(:,1));         % vector of all faces of all cells, true if outer
j = false(6,1);                     % mask, cells can at most have 6 faces,
j(1:4)=true;                        %   extract east, west, north, south
J = j(Gt.cells.faces(:,2));         % vector of faces per cell, true if E,W,N,S
bcIxVE = Gt.cells.faces(I & J, 1);

%{
%% Create figure and plot height
set(0,'CurrentFigure',figure);
plotCellData(G,G.cells.centroids(:,3),'EdgeColor','k','EdgeAlpha',0.05);
set(gca, 'ydir', 'reverse');
view([30 60]), axis tight; drawnow
%}

%% Store data
disp(' -> Writing SleipnerSimple.mat')
if ~isdir(datadir)
   mkdir(datadir);
end
save(fullfile(datadir,'SleipnerSimple'), 'Gt', 'rock', 'rock2D', 'bcIxVE');
end
