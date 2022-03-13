function [G, Gt, rock, rock2D, bcIxVE] = makeSlopingAquiferBig(usemex)
%Make an VE model based upon a data set obtained from the IGEMS project
%
% SYNOPSIS:
%  [G, Gt, bcIx, bcIxVE, rock, rock2D] = makeIGEMSmodel(usemex)
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
% SEE ALSO:
%   `runIGEMS`

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

if nargin<1, usemex=true; end

fname = 'slopingAquiferBig.mat';
try
   disp([' -> Reading ' fname]);
   datadir = fullfile(mrstPath('co2lab'), 'data', 'mat');
   load(fullfile(datadir, fname));
   return;
catch %#ok<*CTCH>
   disp(' -> Reading failed, constructing grid models');
end

%% Read the IGEMS data file
try
   idir = fullfile(mrstPath('co2lab'), 'data', 'igems');
   imodel = 'IGEMS.GRDECL';
   fprintf('    Reading %s\n', imodel);
   grdecl = readGRDECL(fullfile(idir, imodel));
catch
   disp(' -> Download data from: http://www.sintef.no/Projectweb/MRST/');
   disp(['    Putting data in ', idir]);
   untar('https://www.sintef.no/project/MRST/IGEMS.tar.gz',...
      fullfile(mrstPath('co2lab'), 'data'));
   grdecl = readGRDECL(fullfile(idir, imodel));
end


%% Process 3D grid and compute geometry
% First, we map from left-hand to right-hand coordinate system. 
disp(' -> Processing grid');
lines = reshape(grdecl.COORD,6,[]);
lines([2 5],:) = -lines([2 5],:);
grdecl.COORD = lines(:); clear lines

% Then, we process the grid and compute geometry, possibly using
% C-accelerated routines
G = []; 
if usemex,
   mlist = mrstModule;
   mrstModule add libgeometry deckformat;
   try
      G = processgrid(grdecl);
      G = mcomputeGeometry(G);
   catch
      G = [];
   end
   mrstModule('reset',mlist{:});
end

if isempty(G)
   % we either did not want to use mex, or mex failed
   G = processGRDECL(grdecl);
   G = computeGeometry(G);
end

% Adding tags needed by topSurfaceGrid
G.cells.faces = [G.cells.faces, repmat((1:6).', [G.cells.num, 1])];


%% Construct petrophysical model
rock = grdecl2Rock(grdecl);
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


%% Store data
disp([' -> Writing ' fname]);
if ~isdir(datadir)
   mkdir(datadir);
end
save(fullfile(datadir, fname), 'G', 'Gt', 'rock', 'rock2D', 'bcIxVE');
end
