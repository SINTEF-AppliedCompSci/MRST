%% SAIGUP: Solving One-Phase Flow on a Realistic Corner-Point Model
% In the example, we will solve the single-phase pressure equation
%
% $$\nabla\cdot v = q, \qquad v=\textbf{--}\frac{K}{\mu}\nabla p,$$
%
% using the corner-point geometry from synthetic reservoir model from the
% <http://www.nr.no/saigup SAIGUP> study.
%
% The purpose of this example is to demonstrate how the two-point flow solver
% can be applied to compute flow on a real grid model that has degenerate
% cell geometries and non-neighbouring connections arising from a number of
% faults, and inactive cells. A more thorough examination of the model is
% given in a <matlab:edit('showSAIGUP.m') showSAIGUP> example.
%
% Experimentation with this model is continued in
% <matlab:edit('incompNorne2ph.m') incompSAIGUP2ph.m>, in which we solve a
% simplified two-phase model.
mrstModule add incomp

%% Read model data and convert units
% The model data is provided as an ECLIPSE input file that can be read
% using the <matlab:doc('readGRDECL') readGRDECL> function.
grdecl = fullfile(getDatasetPath('SAIGUP'), 'SAIGUP.GRDECL');
grdecl = readGRDECL(grdecl);

%%
% MRST uses the strict SI conventions in all of its internal calculations.
% The SAIGUP model, however, is provided using the ECLIPSE 'METRIC'
% conventions (permeabilities in mD and so on).  We use the functions
% <matlab:doc('getUnitSystem') getUnitSystem> and
% <matlab:doc('convertInputUnits') convertInputUnits> to assist in
% converting the input data to MRST's internal unit conventions.
usys   = getUnitSystem('METRIC');
grdecl = convertInputUnits(grdecl, usys);

%% Define geometry and rock properties
% We generate a space-filling geometry using the
% <matlab:doc('processGRDECL') processGRDECL> function and then compute a
% few geometric primitives (cell volumes, centroids, etc.) by means of the
% <matlab:doc('computeGeometry') computeGeometry> function.
G = processGRDECL(grdecl);
G = computeGeometry(G);

%%
% The media (rock) properties can be extracted by means of the
% <matlab:doc('grdecl2Rock') grdecl2Rock> function.
rock = grdecl2Rock(grdecl, G.cells.indexMap);

%% Modify the permeability to avoid singular tensors
% The input data of the permeability in the SAIGUP realisation is an
% anisotropic tensor with zero vertical permeability in a number of cells.
% We work around this issue by (arbitrarily) assigning the minimum positive
% vertical (cross-layer) permeability to the grid blocks that have zero
% cross-layer permeability.
is_pos                = rock.perm(:, 3) > 0;
rock.perm(~is_pos, 3) = min(rock.perm(is_pos, 3));

%%
% Plot the logarithm of the permeability in the x-direction
clf,
K = log10(rock.perm(:,1));
plotCellData(G,K,'EdgeColor','k','EdgeAlpha',0.1);
set(gca,'dataasp',[15 15 2]);
axis off, view(-110,18);
cs = [0.1 1 10 100 1000 4000];
caxis(log10([min(cs) max(cs)]*milli*darcy));
h=colorbarHist(K,caxis,'South',100); camdolly(.2, 0, 0);
set(h, 'XTick', log10(cs*milli*darcy), 'XTickLabel', num2str(round(cs)'));

%% Set fluid data
% The single fluid has density 1000 kg/m^3 and viscosity 1 cP.
gravity off
fluid      = initSingleFluid('mu' ,    1*centi*poise     , ...
                             'rho', 1000*kilogram/meter^3);

%% Introduce wells
% The reservoir is produced using a set of production wells controlled by
% bottom-hole pressure and rate-controlled injectors. Wells are described
% using a Peacemann model, giving an extra set of equations that need to be
% assembled. For simplicity, all wells are assumed to be vertical and are
% assigned using the logical (i,j) subindex.

% Plot grid outline
clf
plotGrid(G,'FaceColor','none','EdgeAlpha',0.1);
axis tight off, view(-100,8)

% Set six vertical injectors, completed in each layer.
nz = G.cartDims(3);
I = [ 9,  8, 25, 25, 10];
J = [14, 35, 35, 95, 75];
R = [ 2,  4,  4, 1,5]*500*meter^3/day;
radius = .1; W = [];
for i = 1 : numel(I),
   W = verticalWell(W, G, rock, I(i), J(i), 1:nz, 'Type', 'rate', ...
                    'InnerProduct', 'ip_tpf', ...
                    'Val', R(i), 'Radius', radius, 'Comp_i', 1, ...
                    'name', ['I$_{', int2str(i), '}$']);
end
plotGrid(G, vertcat(W.cells), 'FaceColor', 'b','EdgeAlpha',0.1);
prod_off = numel(W);

% Set five vertical producers, completed in each layer.
I = [17, 12, 25, 33, 7];
J = [23, 51, 51, 95, 94];
for i = 1 : numel(I),
   W = verticalWell(W, G, rock, I(i), J(i), 1:nz, 'Type', 'bhp', ...
                    'InnerProduct', 'ip_tpf', ...
                    'Val', 300*barsa(), 'Radius', radius, ...
                    'name', ['P$_{', int2str(i), '}$'], 'Comp_i', 0);
end
plotGrid(G, vertcat(W(prod_off + 1 : end).cells), 'FaceColor', 'r');
plotWell(G,W,'height',30);

%% Compute transmissibilities and solve linear system
% First, we compute one-sided transmissibilities for each local face of
% each cell in the grid. Then, we form a two-point discretization by
% harmonic averaging the one-sided transmissibilities and solve the
% resulting linear system to obtain pressures and fluxes.
T    = computeTrans(G, rock);
rSol = initState(G, W, 350*barsa, 1);
rSol = incompTPFA(rSol, G, T, fluid, 'wells', W);


%%
% Plot the fluxes
clf
cellNo  = rldecode(1:G.cells.num, diff(G.cells.facePos), 2) .';
plotCellData(G, sqrt(accumarray(cellNo,  ...
   abs(convertTo(faceFlux2cellFlux(G, rSol.flux), meter^3/day)))), ...
   'EdgeColor', 'k','EdgeAlpha',0.1);
plotWell(G,W,'height',100,'color','c');
colorbar, axis tight off; view(-80,36)
%%
% Plot the pressure distribution
clf
plotCellData(G, convertTo(rSol.pressure(1:G.cells.num), barsa), ...
   'EdgeColor','k','EdgeAlpha',0.1);
plotWell(G, W, 'height', 100, 'color', 'k');
colorbar, axis tight off; view(-80,36)

%%

% <html>
% <p><font size="-1">
% Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.
% </font></p>
% <p><font size="-1">
% This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).
% </font></p>
% <p><font size="-1">
% MRST is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% </font></p>
% <p><font size="-1">
% MRST is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% </font></p>
% <p><font size="-1">
% You should have received a copy of the GNU General Public License
% along with MRST.  If not, see
% <a href="http://www.gnu.org/licenses/">http://www.gnu.org/licenses</a>.
% </font></p>
% </html>
