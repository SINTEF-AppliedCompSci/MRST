%% Simulate a transport problem with the MsRSB method
% This example is a simple demonstration in how to use the multiscale
% solver to solve a transport problem using the incomp module.

mrstModule add mrst-gui spe10 coarsegrid incomp


%%  Set up the fine-scale model
% We set up a layer of SPE10, just as in the introductory single-phase
% <matlab:edit('exampleMs.m') exampleMs> example.

% Grid an petrophysics
layerNo = 85;
[G, ~, rock] = getSPE10setup(layerNo);
mp = 0.01;
rock.poro(rock.poro < mp) = mp;

% Transmissibilities on fine scale
T = computeTrans(G, rock);

% Inject 1 PV over 10 years
time = 10*year;
Nt = 100;
pv = poreVolume(G, rock);

% Set up wells in each corner
injr = sum(pv)/(time);
W = [];
W = verticalWell(W, G, rock, 1, 1, [], 'type', 'rate', ...
    'val', injr, 'comp_i', [1, 0]);
W = verticalWell(W, G, rock, G.cartDims(1), G.cartDims(2), [],...
    'type', 'bhp' , 'val', 100*barsa, 'comp_i', [1, 0]);

% Plot the permeability field (log-scale).
figure(1); clf
plotToolbar(G, log10(rock.perm(:, 1)));
axis equal tight off
colormap jet
view(90, 90);

%% Define coarse partition and support regions
% Set up a uniform coarsegrid
mrstModule add msrsb

coarsen = [5 10 5];
coarsedims = ceil(G.cartDims./coarsen);

% Generate partition vector
p = partitionUI(G, coarsedims);
% Grid structure from partition vector
CG = generateCoarseGrid(G, p);
% Add centroids / geometry information on coarse grid
CG = coarsenGeometry(CG);
% Store the support of each cell (used by multiscale basis construction)
CG = storeInteractionRegionCart(CG);

%% Set up basis functions
A = getIncomp1PhMatrix(G, T);
basis = getMultiscaleBasis(CG, A, 'type', 'msrsb');

%% Simulate the problem
% We simulate the whole injection period with a fluid model with equal
% viscosities and quadratic Corey coefficients for the relative
% permeability. Since the total mobility does not change significantly, we
% opt to keep the basis functions static throughout the simulation.
% Moreover, we only use a single multiscale solve and do not introduce
% extra iterations to reduce the fine-scale residual towared machine
% precision. As a result, the multiscale approximation will deviate
% somewhat from the fine-scale TPFA solution.

dt = time/Nt;
figure(1);
set(gcf, 'position', [50, 50, 1000, 500])

fluid = initSimpleFluid('mu', [1, 1]*centi*poise, 'n', [2, 2], 'rho', [0, 0]);

state0 = initResSol(G, 0, [0, 1]);
gravity reset off

psolve = @(state) incompTPFA(state, G, T, fluid, 'Wells', W);
mssolve = @(state) incompMultiscale(state, CG, T, fluid, basis, 'Wells', W);
solver = @(state) implicitTransport(state, G, dt, rock, fluid, 'wells', W);

states = psolve(state0);
states_ms = mssolve(state0);

h = waitbar(0, 'Starting simulation...');
for i = 1:Nt
    waitbar(i/Nt, h, sprintf('Step %d of %d: ', i, Nt))
    
    state = psolve(states(end));
    states = [states; solver(state)];
    
    state = mssolve(states_ms(end));
    states_ms = [states_ms; solver(state)];
    
    s_ms = states_ms(end).s(:, 1);
    s_ref = states(end).s(:, 1);
    figure(1), clf
    subplot(1, 3, 1)
    plotCellData(G, s_ref, 'EdgeColor', 'none')
    axis tight off
    title('Reference');
    
    subplot(1, 3, 2)
    plotCellData(G, s_ms, 'EdgeColor', 'none')
    axis tight off
    title('MsRSB')
    
    subplot(1, 3, 3)
    plotCellData(G, abs(s_ref - s_ms), 'EdgeColor', 'none')
    axis tight off
    caxis([0, 1])
    title('Saturation error')
    drawnow
end
close(h)

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
