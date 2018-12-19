%% Use of Peacemann well models
% In this example we will demonstrate how to set up a flow problems with
% two wells, one rate-controlled, vertical well and one horizontal well
% controlled by bottom-hole pressure. The reservoir is a regular box with
% homogeneous petrophysical properties.

mrstModule add incomp

%% Set up reservoir model
[nx,ny,nz] = deal(20,20,5);
G = computeGeometry( cartGrid([nx,ny,nz], [500 500 25]) );
rock  = makeRock(G, 100.*milli*darcy, .2);
fluid = initSingleFluid('mu', 1*centi*poise,'rho', 1014*kilogram/meter^3);
hT    = computeTrans(G, rock);

%% Add wells wells
W = verticalWell([], G, rock, 1, 1, 1:nz, 'Type', 'rate', 'Comp_i', 1,...
                'Val', 3e3/day, 'Radius', .12*meter, 'name', 'I');
disp('Well #1: '); display(W(1));

W = addWell(W, G, rock, nx : ny : nx*ny, 'Type', 'bhp', 'Comp_i', 1, ...
            'Val', 1.0e5, 'Radius', .12*meter, 'Dir', 'y', 'name', 'P');
disp('Well #2: '); display(W(2));

state = initState(G, W, 0);

%%
% We plot the wells to check if the wells are placed as we wanted them.
% (The plot will later be moved to subplot(2,2,1), hence we first find the
% corresponding axes position before generating the handle graphics).
subplot(2,2,1), pos = get(gca,'Position'); clf
plotGrid(G, 'FaceColor', 'none');
view(3), camproj perspective, axis tight off,
plotWell(G, W(1), 'radius',  1, 'height', 5, 'color', 'r');
plotWell(G, W(2), 'radius', .5, 'height', 5, 'color', 'b');

%% Assemble and solve system
gravity reset on;
state = incompTPFA(state, G, hT, fluid, 'wells', W);

%% Report results
% We move the plot of the grids and wells to the upper-left subplot. The
% producer inflow profile is shown in the upper-right and the cell
% pressures in the lower-left subplot. In the lower-right subplot, we show
% the flux intensity, which must be constructed by averaging over cell
% faces

%subplot(2,2,1)
   set(gca, 'Position', pos);  % move the current plot

subplot(2,2,2)
   plot(convertTo(-state.wellSol(2).flux, meter^3/day),'o')
   title('Producer inflow profile [m^3/d]');

subplot(2,2,3)
   plotCellData(G, convertTo(state.pressure(1:G.cells.num), barsa),'EdgeAlpha',.1);
   title('Pressure [bar]')
   view(3), camproj perspective, axis tight off

subplot(2,2,4)
   [i j k] = ind2sub(G.cartDims, 1:G.cells.num);
   I = false(nx,1); I([1 end])=true;
   J = false(ny,1); J(end)=true;
   K = false(nz,1); K([1 end]) = true;
   cf = accumarray(getCellNoFaces(G), ...
      abs(faceFlux2cellFlux(G, state.flux)));
   plotCellData(G, convertTo(cf, meter^3/day), I(i) | J(j) | K(k),'EdgeAlpha',.1);
   title('Flux intensity [m^3/day]')
   view(-40,20), camproj perspective, axis tight, box on
   set(gca,'XTick',[],'YTick',[],'ZTick',[]);

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
