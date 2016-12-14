function [G, Gt, rock, rock2D, bcIxVE] = makeSleipnerVEmodel(varargin)
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
%
% OPTIONS:
%   usemex          - Flag: if true, use C-accelerated routines for
%                   processing Eclipse input and computing geometry
%   assign_coords   - Flag: if true, will read physical coordinates from
%                   Sleipner's M9X1.grdecl datafile, and save as 
%                   SleipnerGlobalCoords.mat, otherwise Sleipner.mat
%
% DESCRIPTION:
%
% SEE ALSO:
%   runSleipner

opt.usemex          = true;
opt.assign_coords   = false;
opt = merge_options(opt, varargin{:});

try
   datadir = fullfile(mrstPath('co2lab'), 'data', 'mat');
   if ~opt.assign_coords
       disp(' -> Reading Sleipner.mat');
       load(fullfile(datadir,'Sleipner'));
       return;
   else
       disp(' -> Reading SleipnerGlobalCoords.mat');
       load(fullfile(datadir,'SleipnerGlobalCoords'));
       return;
   end
catch %#ok<*CTCH>
   disp('    Data set has not yet been constructed.');
end

%% Read data
% First loading: to get cartDims, COORD, ZCORN, ACTNUM, PERMX, PERMZ, PORO
try
   sdir = fullfile('data', 'sleipner');
   disp([' -> Reading data from: ' sdir]);
   grdecl = readGRDECL(fullfile(mrstPath('co2lab'), sdir, 'SLEIPNER.DATA'));
catch
   fprintf(1, '    Reading of SLEIPNER.DATA failed, please dowload data manually');
   fprintf(1, ' following instructions\n    in "%s"\n', fullfile(sdir,'README'));
   return;
end

% Second loading (optional): to get MAPUNITS, MAPAXES
if opt.assign_coords
   try
       moduleCheck('deckformat', 'mex');
       sl_file = fullfile(mrstPath('co2lab'), 'data', 'sleipner', 'M9X1.grdecl'); % IEAGHG 
       fn      = fopen(sl_file);
       gr      = readGRID(fn, fileparts(sl_file), initializeDeck());
       fclose(fn); 
       % Add MAPAXES and MAPUNITS to grdecl
       grdecl.MAPAXES = gr.GRID.MAPAXES;
       grdecl.MAPUNITS = gr.GRID.MAPUNITS;
       clear gr sl_file
   catch
       fprintf(1, '    Reading of M9X1.grdecl failed, please dowload data manually');
       fprintf(1, ' following instructions\n    in "%s"\n', fullfile(sdir,'README'));
       return;
   end
end

%% Process 3D grid and compute geometry 
disp(' -> Processing grid');

if ~opt.assign_coords
    % First, we map from left-hand to right-hand coordinate system. 
    lines = reshape(grdecl.COORD,6,[]);
    lines([2 5],:) = -lines([2 5],:);
    grdecl.COORD = lines(:); clear lines
else
    % First, we recompute X and Y coordinates in terms of the provided axes
    % (depths, Z, do not require any recomputation)
    coords        = reshape(grdecl.COORD,3,[])';
    coords(:,1:2) = mapAxes(coords(:,1:2), grdecl.MAPAXES);
    coords        = coords';
    grdecl.COORD  = coords(:); clear coords
end

% Then, we remove the bottom and top layers that contain shale
grdecl.ACTNUM(grdecl.PERMX<200) = 0;

% Next, we process the grid and compute geometry, possibly using
% C-accelerated routines
G = [];
if opt.usemex,
   mlist = mrstModule;
   mrstModule add libgeometry opm_gridprocessing
   try
      G = processgrid(grdecl);
      G = mcomputeGeometry(G);
   catch
      G = [];
   end
   mrstModule('reset', mlist{:});
end
if isempty(G)
   % We either did not want to use mex, or mex failed
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
if ~isdir(datadir)
    mkdir(datadir);
end
if ~opt.assign_coords
    disp(' -> Writing Sleipner.mat')
    save(fullfile(datadir,'Sleipner'), 'G', 'Gt', 'rock', 'rock2D', 'bcIxVE');
else
    disp(' -> Writing SleipnerGlobalCoords.mat')
    save(fullfile(datadir,'SleipnerGlobalCoords'), 'G', 'Gt', 'rock', 'rock2D', 'bcIxVE');
end
end
