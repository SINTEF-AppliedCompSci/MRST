%% Example 2: Grid-orientation effect
% The multipoint flux-approximation (MPFA-O) scheme is developed to be
% consistent on grids that are not necessarily K-orthogonal. This example
% demonstrates the basic use of the MPFA-O pressure solver by applying it
% to a single-phase flow problem posed on a square [0,30]x[0,30] m^2 with
% homogeneous and isotropic permeability. To discretize the problem, we
% introduce a skewed, curvilinear grid in which the majority of the grid
% cells are not K-orthogonal. The classical TPFA scheme can therefore be
% expected to give significant grid-orientation effects. To investigate
% this, we consider a single-phase flow problem with a prescribed pressure
% drop from the left to the right boundary, which gives an analytical
% solution p(x,y)=2-x/30. We compare the MPFA-O solution with solutions
% computed by the TPFA method and a mimetic method.

mrstModule add incomp mimetic mpfa

%% Set up simulation model
gravity off
G = cartGrid([30, 30]);
G.nodes.coords = twister(G.nodes.coords);
G = computeGeometry(G);

rock.perm = 0.1*darcy*ones(G.cells.num, 1);

bc = pside([], G, 'left',  2);
bc = pside(bc, G, 'right', 1);

fluid = initSingleFluid('mu' ,    1*centi*poise     , ...
                        'rho', 1014*kilogram/meter^3);

state0 = initResSol(G, 0, 0);
%% MPFA-O method
fprintf('MPFA-O method\t... ')
tic
T1  = computeMultiPointTrans(G, rock);
xr1 = incompMPFA(state0, G, T1, fluid, 'bc', bc, 'MatrixOutput', true);
toc

%% Mimetic method
fprintf('Mimetic method\t... ')
tic
S   = computeMimeticIP(G, rock);
xr2 = incompMimetic(state0, G, S, fluid, 'bc', bc);
toc

%% TPFA method
fprintf('TPFA Method\t... ')
tic
T2  = computeTrans(G, rock);
xr3 = incompTPFA(state0, G, T2, fluid, 'bc', bc,'MatrixOutput', true);
toc

%% Plot solutions
hf2cn      = gridCellNo(G);
flux_int   = @(x) accumarray(hf2cn, abs(x.flux(G.cells.faces(:,1))));
plot_var   = @(v) plotCellData(G, v, 'EdgeColor','none');
plot_press = @(x) plot_var(x.pressure(1:G.cells.num));
plot_flux  = @(x) plot_var(convertTo(flux_int(x), meter^3/day));

clf, set(gcf,'Position',[300 250 1000 500]);

subplot(2,3,1),
plot_flux(xr1); cax=caxis();  axis equal tight, title('Flux: MPFA-O')

subplot(2,3,2),
plot_flux(xr2); caxis(cax); axis equal tight, title('Flux: Mimetic')

subplot(2,3,3),
plot_flux(xr3); caxis(cax);  axis equal tight, title('Flux: TPFA')
colorbar('Position',[.92 .58 .02 .34])

subplot(2,3,4),
plot_press(xr1); cax=caxis();  axis equal tight, title('Pressure: MPFA-O')

subplot(2,3,5),
plot_press(xr2); caxis(cax); axis equal tight, title('Pressure: Mimetic')

subplot(2,3,6),
plot_press(xr3); caxis(cax);  axis equal tight, title('Pressure: TPFA')
colorbar('Position',[.92 .11 .02 .34])

%% Compute discrepancies in flux and errors in pressure
p = struct();
p.pressure = 2 - G.cells.centroids(:,1)/G.cartDims(1);
err_press  = @(x1, x2) ...
    norm(x1.pressure - x2.pressure, inf) / norm(x1.pressure, inf);
err_flux   = @(x1, x2) norm(flux_int(x1) - flux_int(x2), inf);

fprintf(['\nMaximum difference in face fluxes:\n', ...
         '\to MPFA-O /TPFA   : %.15e\n',   ...
         '\to MPFA-O /Mimetic: %.15e\n',   ...
         '\to Mimetic/TPFA   : %.15e\n\n', ...
         ], ...
        err_flux(xr1, xr3), err_flux(xr1, xr2), err_flux(xr2, xr3));

fprintf(['Relative error in cell pressures:\n', ...
         '\to MPFA-O         : %.15e\n',  ...
         '\to Mimetic        : %.15e\n',  ...
         '\to TPFA           : %.15e\n',  ...
         ], ...
         err_press(xr1, p), err_press(xr2, p), err_press(xr3, p));

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
