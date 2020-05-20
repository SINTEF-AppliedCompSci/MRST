%% Pressure Solver: Example of a realistic Field Model
% In the example, we will solve the single-phase, incompressible pressure
% equation using the corner-point geometry from synthetic reservoir model
% from the SAIGUP study.
%
% The purpose of this example is to demonstrate how the two-point flow solver
% can be applied to compute flow on a real grid model that has degenerate
% cell geometries and non-neighbouring connections arising from a number of
% faults, and inactive cells.

%% Check for existence of input model data
% The model can be downloaded from the the MRST page
%
% http://www.sintef.no/Projectweb/MRST/
grdecl = fullfile(getDatasetPath('SAIGUP'), 'SAIGUP.GRDECL');
if ~exist(grdecl, 'file'),
   error('SAIGUP model data is not available.')
end

%% Load and process grid model
% The model data is provided as an ECLIPSE input file. MRST uses the strict
% SI conventions in all of its internal calculations. The SAIGUP model,
% however, is provided using the ECLIPSE 'METRIC' conventions
% (permeabilities in mD and so on) and we must therefore convert the input
% data to MRST's internal unit conventions.
grdecl = readGRDECL(grdecl);
usys   = getUnitSystem('METRIC');
grdecl = convertInputUnits(grdecl, usys);
G = processGRDECL(grdecl);
G = computeGeometry(G);

% To speed up the processing, one can use a C-accelerated versions of the
% gridprocessing routines (provided you have a compatible compiler):
% mrstModule add libgeometry deckformat
% G = mcomputeGeometry(processgrid(grdecl));

%% Get petrophysical data
% The input data of the permeability in the SAIGUP realisation is an
% anisotropic tensor with zero vertical permeability in a number of cells.
% As a result some parts of the reservoir to be completely sealed from
% the wells. This will cause problems for the time-of-flight solver, which
% requires that all cells in the model must be flooded after some finite
% time that can be arbitrarily large. We work around this issue by
% assigning a small constant times the minimum positive vertical
% (cross-layer) permeability to the grid blocks that have zero cross-layer
% permeability.
rock = grdecl2Rock(grdecl, G.cells.indexMap);
is_pos = rock.perm(:, 3) > 0;
rock.perm(~is_pos, 3) = 1e-6*min(rock.perm(is_pos, 3));

hT   = computeTrans(G, rock);

%% Set fluid data
mrstModule add incomp
gravity reset on
fluid = initSingleFluid('mu',1*centi*poise,'rho', 1000*kilogram/meter^3);

%% Introduce wells
% The reservoir is produced using a set of production wells controlled by
% bottom-hole pressure and rate-controlled injectors. Wells are described
% using a Peacemann model, giving an extra set of equations that need to be
% assembled. For simplicity, all wells are assumed to be vertical and are
% assigned using the logical (i,j) subindex.

% Plot grid outline
figure('position',[440 317 866 480]);
plotCellData(G,log10(rock.perm(:,1)), ...
   'EdgeColor','k','EdgeAlpha',.1,'FaceAlpha',.5);
axis tight off, view(-100,20)

% Set eight vertical injectors around the perimeter of the model, completed
% in each layer.
nz = G.cartDims(3);
I = [ 3, 20,  3, 25,  3, 30,  5, 29];
J = [ 4,  3, 35, 35, 70, 70,113,113];
R = [ 1,  3,  3,  3,  2,  4,  2,  3]*500*meter^3/day;
W = [];
for i = 1 : numel(I),
   W = verticalWell(W, G, rock, I(i), J(i), 1:nz, 'Type', 'rate', ...
      'Val', R(i), 'Radius', .1*meter, 'Comp_i', 1, ...
      'name', ['I$_{', int2str(i), '}$']);
end
plotWell(G, W, 'height', 30, 'color', 'k');
in = numel(W);

% Set six vertical producers, completed in each layer.
I = [15, 12, 25, 21, 29, 12];
J = [25, 51, 51, 60, 95, 90];
for i = 1 : numel(I),
   W = verticalWell(W, G, rock, I(i), J(i), 1:nz, 'Type', 'bhp', ...
      'Val', 200*barsa(), 'Radius', .1*meter, ...
      'name', ['P$_{', int2str(i), '}$'], 'Comp_i',1);
end
plotWell(G,W(in+1:end),'height',30,'color','b');

%% Assemble and solve system, plot results
state = initState(G, W, 350*barsa, 1);
state = incompTPFA(state, G, hT, fluid, 'wells', W);

figure('position',[440 317 866 480]);
plotCellData(G, convertTo(state.pressure(1:G.cells.num), barsa), ...
   'EdgeColor','k','EdgeAlpha',0.1);
plotWell(G, W(1:in),     'height', 100, 'color', 'b');
plotWell(G, W(in+1:end), 'height', 100, 'color', 'k');
axis off; view(-80,36)
h=colorbar; set(h,'Position',[0.88 0.15 0.03 0.67]);

%% Time-of-flight analysis
mrstModule add diagnostics
tf = computeTimeOfFlight(state, G, rock, 'wells', W)/year;
tb = computeTimeOfFlight(state, G, rock, 'wells', W, 'reverse', true)/year;

figure('position',[440 317 866 480]);
plotCellData(G,tf+tb,tf+tb<50,'EdgeColor','k','EdgeAlpha',0.1);
plotWell(G, W(1:in),     'height', 100, 'color', 'b');
plotWell(G, W(in+1:end), 'height', 100, 'color', 'k');
plotGrid(G,'FaceColor','none','edgealpha',.05);
axis off; view(-80,36)
caxis([0 50]);
h=colorbar; set(h,'Position',[0.88 0.15 0.03 0.67]);

%{
Copyright 2009-2018 SINTEF Digital, Mathematics & Cybernetics.

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
