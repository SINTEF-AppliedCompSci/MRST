function [G, Gt, rock, rock2D, bcIxVE] = makeSleipnerVEmodel(varargin)
%Make a VE model based upon the Sleipner data set from ieaghg.org
%
% SYNOPSIS:
%   [G, Gt, rock, rock2D, bcIxVE] = makeSleipnerVEmodel()
%   [G, Gt, rock, rock2D, bcIxVE] = makeSleipnerVEmodel('pn1', pv1, ...)
%
% PARAMETERS:
%   'pn'/pv - List of 'key'/value pairs defining optional parameters.  The
%             supported options are:
%
%               usemex --
%                   Whether or not to use C-accelerated routines for
%                   processing ECLIPSE input and computing geometry.
%
%               assign_coords --
%                   Whether or not to read physical coordinates from
%                   Sleipner's M9X1.grdecl datafile.  If set, do read
%                   coordinates and save as 'SleipnerGlobalCoords.mat',
%                   otherwise save result as 'Sleipner.mat'.
%
% RETURNS:
%   G      - Data structure for 3D grid
%   Gt     - Data structure for topsurface grid
%   rock   - Data structure for 3D rock parameters
%   rock2D - Data structure for rock parameters for topsurface grid
%   bcIxVE - Index for pressure boundary conditions in topsurface grid
%
% SEE ALSO:
%   `runSleipner`.

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}

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

% Read data
% First loading: to get cartDims, COORD, ZCORN, ACTNUM, PERMX, PERMZ, PORO
sdir = fullfile(mrstDataDirectory, 'Sleipner');
try
   disp([' -> Reading data from: ' sdir]);
   grdecl = readGRDECL(fullfile(sdir, 'SLEIPNER.DATA'));
catch
   error('Dataset:Missing', ...
        ['Reading of SLEIPNER.DATA failed, please download data ', ...
         'manually following instructions in\n    ''%s'''], ...
         fullfile(sdir,'README'));
end

% Second loading (optional): to get MAPUNITS, MAPAXES
if opt.assign_coords
   try
       moduleCheck('deckformat', 'mex');
       sl_file = fullfile(sdir, 'M9X1.grdecl'); % IEAGHG
       fn      = fopen(sl_file);
       gr      = readGRID(fn, fileparts(sl_file), initializeDeck());
       fclose(fn);
       % Add MAPAXES and MAPUNITS to grdecl
       grdecl.MAPAXES = gr.GRID.MAPAXES;
       grdecl.MAPUNITS = gr.GRID.MAPUNITS;
       clear gr sl_file
   catch
       error('Dataset:Missing', ...
            ['Reading of M9X1.grdecl failed. Data is missing from %s.'], ...
            sdir);
   end
end

% Process 3D grid and compute geometry
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
   mrstModule add libgeometry deckformat
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


% Construct petrophysical model
rock = grdecl2Rock(grdecl, G.cells.indexMap);
rock.perm = convertFrom(rock.perm, milli*darcy);
clear grdecl


% Construct top-surface grid
disp(' -> Constructing top-surface grid');
[Gt, G] = topSurfaceGrid(G);
rock2D  = averageRock(rock, Gt);


% Find pressure boundary
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

% Store data
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
