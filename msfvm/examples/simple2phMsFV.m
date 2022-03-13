%% A simple two phase problem solved using the Multiscale Finite Volume method
% The multiscale finite volume method can easily be used in place of a
% regular pressure solver for incompressible transport. This example
% demonstrates a two phase solver on a 2D grid.

mrstModule add coarsegrid msfvm incomp msrsb
%% Construct simple 2D Cartesian test case
nx = 50; ny = 50;
Nx = 5; Ny = 5;
G         = cartGrid([nx, ny], [1000, 1000]);
G         = computeGeometry(G);

% Plot each timestep
doPlot = true;

p  = partitionUI(G, [Nx, Ny]);
CG = generateCoarseGrid(G, p);
CG = coarsenGeometry(CG);
%% Generate dual grid
DG = partitionUIdual(CG, [Nx, Ny]);
CG.dual = DG;
%% Visualize
clf;
plotDual(G, DG)
%% Uniform permeability
rock = makeRock(G, 100*milli*darcy, 0.3);
T = computeTrans(G, rock);

%% Define a simple 2 phase fluid
fluid = initSimpleFluid('mu' , [   1,  10]*centi*poise     , ...
                        'rho', [1014, 859]*kilogram/meter^3, ...
                        'n'  , [   2,   2]);


%% Setup a producer / injector pair of wells
rate = 0.01*meter^3/second;
bhp  = 100*barsa;
radius = 0.05;
% Injector in lower left corner
W = [];
W = verticalWell(W, G, rock, 5, 5, [],      ...
            'Type', 'rate' , 'Val', rate, ...
            'Radius', radius, 'InnerProduct', 'ip_tpf', ...
            'Comp_i', [1, 0]);
% Producer in upper right corner
W = verticalWell(W, G, rock, nx - 5, ny - 5, [],     ...
            'Type', 'bhp' , 'Val', bhp, ...
            'Radius', radius, 'InnerProduct', 'ip_tpf', ...
            'Comp_i', [0, 1]);
%% Set up solution structures with only one phase
refSol    = initState(G, W, 0, [0, 1]);
msSol     = initState(G, W, 0, [0, 1]);

gravity off
verbose = false;

%% Set up pressure and transport solvers
% Reference TPFA
r_psolve = @(state) incompTPFA(state, G, T, fluid, 'wells', W);
% MsFV solver. We use the more modern version of the solver found in the
% msrsb module. From the linear system A, we construct a set of msfv basis
% functions.
A = getIncomp1PhMatrix(G, T, msSol, fluid);
basis = getMultiscaleBasis(CG, A, 'type', 'msfv');
psolve = @(state) incompMultiscale(state, CG, T, fluid, basis, 'Wells', W);
% Implicit transport solver
tsolve   = @(state, dT) implicitTransport(state, G, dT, rock, ...
                                                fluid, 'wells', W, ...
                                                'verbose', verbose);
%%
% Alternatively we could have defined an explicit transport solver by
%
% tsolve = @(state, dT, fluid) explicitTransport(state, G, dT, rock, fluid, ...
%                                        'wells', W, 'verbose', verbose);

%% Solve initial pressure in reservoir
% We solve and plot the pressure in the reservoir at t=0.
refSol= r_psolve(refSol);
msSol = psolve(msSol);
subplot(2,1,1)
plotCellData(G, refSol.pressure); axis tight; colorbar;
title('Pressure Ref')
cbar = caxis();
subplot(2,1,2)
plotCellData(G, msSol.pressure); axis tight; colorbar;
title('Pressure MS')
caxis(cbar)

%% Transport loop
% We solve the two-phase system using a sequential splitting in which the
% pressure and fluxes are computed by solving the flow equation and then
% held fixed as the saturation is advanced according to the transport
% equation.
time   = 1*year;
dT     = time/60;

%% Start the main loop
% Iterate through the time steps and plot the saturation profiles along the
% way.
t = 0;
while t < time
    % Solve transport equations using the transport solver
    msSol  = tsolve(msSol , dT);
    refSol = tsolve(refSol, dT);

    % Update the pressure based on the new saturations
    msSol    = psolve(msSol);
    refSol   = r_psolve(refSol);
    % Increase time and continue if we do not want to plot saturations
    if doPlot
        clf;
        % Saturation plot
        subplot(2,2,1)
        plotGrid(G, 'FaceColor', 'None', 'EdgeAlpha', 0)
        plotCellData(G, refSol.s(:,1), refSol.s(:,1) > 1e-4); axis tight; colorbar;
        title('Saturation Ref')
        caxis([0 1]);
        subplot(2,2,2)
        plotGrid(G, 'FaceColor', 'None', 'EdgeAlpha', 0)
        plotCellData(G, msSol.s(:,1), msSol.s(:,1) > 1e-4); axis tight; colorbar;
        title('Saturation MSFV')
        % Align colorbars
        caxis([0 1])
        % Pressure plot
        subplot(2,2,3)
        plotCellData(G, refSol.pressure); axis tight; colorbar;
        title('Pressure Ref')
        minP = min(min(refSol.pressure), min(msSol.pressure));
        maxP = max(max(refSol.pressure), max(msSol.pressure));
        cbar = [minP, maxP];
        subplot(2,2,4)
        hs = plotCellData(G, msSol.pressure); axis tight; colorbar;
        title('Pressure MSFV')
        caxis(cbar)
        drawnow
    end
    t = t + dT;
end

%%
% <html>
% <p><font size="-1">
% Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.
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
