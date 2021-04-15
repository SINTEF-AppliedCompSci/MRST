function [G, rock, grdecl, Gfull, rockfull] = getMultiLayerSleipnerGrid(varargin)
%Get new multilayered benchmark Sleipner grid
%
% SYNOPSIS:
%   [G, rock, grdecl, Gfull, rockfull] = getMultiLayerSleipnerGrid()
%   [G, rock, grdecl, Gfull, rockfull] = getMultiLayerSleipnerGrid('coarseDimXY',[nx ny])
%
%
% RETURNS:
%   G         - grid with caprock cells removed.
%   rock      - rock corresponding to G. 
%   grdecl    - mrst grdecl structure.
%   Gfull     - grid with caprock cells.
%   rockfull  - rock corresponding to Gfull.
%
%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

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

opt = struct('coarseDimXY',[]);
opt = merge_options(opt, varargin{:});

datafolder = fullfile(getDatasetPath('sleipner2019'),...
    'Sleipner_Reference_Model_2019_Grid','data');

grdecl = readGRDECL(fullfile(datafolder, 'Sleipner_Reference_Model.grdecl'));

% Reshaping
lines = reshape(grdecl.COORD,6,[]);
grdecl.COORD = lines(:); clear lines

% Add PERMY and PERMZ
grdecl.PERMY = grdecl.PERMX;
grdecl.PERMZ = grdecl.PERMX;


% 
if ~isempty(opt.coarseDimXY)
    require mrst-experimental
    coarseDim = [opt.coarseDimXY 1];
    grdecl = cutGrdecl(grdecl,[1 64; 1 112; 1 263]);
    grdecl = coarseGrdecl(grdecl, coarseDim,'only_grid',false);
end

% Second loading of Sleipner Eclispe grid, to get MAPAXES
require deckformat
sl_file = fullfile(datafolder, 'Sleipner_Reference_Model.grdecl');
fn      = fopen(sl_file);
gr  = readGRID(fn, fileparts(sl_file), initializeDeck());
% this grdecl contains: GRID, and others. grdecl.GRID contains
% MAPUNITS, MAPAXES, cartDims, COORD, ZCORN, ACTNUM
fclose(fn);


% Add data loaded from first loading of Sleipner Eclispe grid
grdecl.MAPAXES = gr.GRID.MAPAXES;
grdecl.MAPUNITS = gr.GRID.MAPUNITS;
clear gr sl_file

%% Transform coordinates

% If required, recompute X and Y coordinates in terms of
% the provided axes (depths, Z, do not require any
% recomputation)
coords        = reshape(grdecl.COORD,3,[])';
coords(:,1:2) = mapAxes(coords(:,1:2), grdecl.MAPAXES);
coords        = coords';
grdecl.COORD  = coords(:); clear coords


%% Next, we process the grid and compute geometry
require libgeometry
Gfull = processGRDECL(grdecl); % processgrid() didn't work
Gfull = mcomputeGeometry(Gfull);

rockfull        = grdecl2Rock(grdecl, Gfull.cells.indexMap);
rockfull.perm   = convertFrom(rockfull.perm, milli*darcy);

%% Remove caprock cells.
% First 10 layers of the model correspond to the caprock. (See
% Final_Report_04012019.pdf Fig. 42)

caprockCells = Gfull.cartDims(1)*Gfull.cartDims(2)*10;

[G,gc] = extractSubgrid(Gfull,caprockCells+1:1:Gfull.cells.num);

rock.perm = rockfull.perm(gc,:);
rock.poro = rockfull.poro(gc,:);







