%% Compare grid-orientation effects for TPFA/mimetic schemes

mrstModule add incomp mimetic streamlines diagnostics

% Rectangular reservoir with a skew grid.
G = cartGrid([41,20],[2,1]);
makeSkew = @(c) c(:,1) + .4*(1-(c(:,1)-1).^2).*(1-c(:,2));
G.nodes.coords(:,1) = 2*makeSkew(G.nodes.coords);
% G.nodes.coords = twister(G.nodes.coords);
% G.nodes.coords(:,1) = 2*G.nodes.coords(:,1);
G = computeGeometry(G);

% Homogeneous reservoir properties
rock = makeRock(G, 100*milli*darcy, .2);
pv   = sum(poreVolume(G,rock));

% Symmetric well pattern
srcCells = findEnclosingCell(G,[2 .975; .5 .025; 3.5 .025]);
src = addSource([], srcCells, [pv; -.5*pv; -.5*pv]);

% Single-phase fluid
fluid = initSingleFluid('mu', 1*centi*poise,'rho', 1000*kilogram/meter^3);

%%
% Figure and figure settings
figure('Position',[400 460 900 350]);
parg = {'EdgeColor','k','EdgeAlpha',.05};

%%
% TPFA solution
hT   = computeTrans(G, rock);
s_tp = initState(G,[], 0);
s_tp = incompTPFA(s_tp, G, hT, fluid, 'src', src);
subplot(2,2,1);
plotCellData(G, s_tp.pressure, 'EdgeColor', 'k', 'EdgeAlpha', .05); 

%%
% Mimetic solution
S = computeMimeticIP(G, rock);
s_mi = initState(G, [], 0);
s_mi = incompMimetic(s_mi, G, S, fluid, 'src', src);
subplot(2,2,2);
plotCellData(G, s_mi.pressure, 'EdgeColor', 'k', 'EdgeAlpha', .05); 

%%
% Compare time-of-flight
subplot(2,2,3);
tof_tp = computeTimeOfFlight(s_tp, G, rock, 'src', src);
plotCellData(G, tof_tp, tof_tp<.2,'EdgeColor','none'); caxis([0 .2]); box on
seed = floor(G.cells.num/5)+(1:G.cartDims(1))';
hf = streamline(pollock(G, s_tp, seed, 'substeps', 1) );
hb = streamline(pollock(G, s_tp, seed, 'substeps', 1, 'reverse' , true));
set ([ hf ; hb ], 'Color' , 'k' );

subplot(2,2,4);
tof_mi = computeTimeOfFlight(s_mi, G, rock, 'src', src);
plotCellData(G, tof_mi, tof_mi<.2,'EdgeColor','none'); caxis([0 .2]); box on
hf = streamline(pollock(G, s_mi, seed, 'substeps', 1) );
hb = streamline(pollock(G, s_mi, seed, 'substeps', 1, 'reverse' , true));
set ([ hf ; hb ], 'Color' , 'k' );

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
