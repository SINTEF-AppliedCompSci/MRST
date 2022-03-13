%% Water injection in a field model using the MsRSB-method
% This example demonstrates the use of the MsRSB method to a synthethic
% waterflood on the grid and petrophysical properties from a real field. 
%
% This example is a modified version of Example 4.3.6 in
%  A multiscale restriction-smoothed basis method for high contrast porous
%  media represented on unstructured grids. J. Comput. Phys, Vol. 304, pp.
%  46-71, 2016. DOI: 10.1016/j.jcp.2015.10.010
%
% Note that this example uses MEX-accelerated computation of basis functions
% which requires that MEX is set up and configured to work with an
% installed C++ compiler on your machine.

mrstModule add coarsegrid msrsb ad-core mrst-gui incomp

%% Read and set up the Norne model
% We use a subset of the Norne model that is used in other examples in
% MRST. We set up a model with anisotropic permeability based on the
% vertical permeability.
mrstModule add deckformat

if ~(makeNorneSubsetAvailable() && makeNorneGRDECL())
   error('Unable to obtain simulation model subset');
end

grdecl = fullfile(getDatasetPath('norne'), 'NORNE.GRDECL');
grdecl = readGRDECL(grdecl);
usys   = getUnitSystem('METRIC');
grdecl = convertInputUnits(grdecl, usys);

G = processGRDECL(grdecl, 'PreserveCpNodes', true);
G = computeGeometry(G(1));
rock = grdecl2Rock(grdecl, G.cells.indexMap);

% Backwards compatible plotting
if isnumeric(gcf)
    myZoom = @zoom;
else
    myZoom = @(varargin) [];
end
%% Plot permeability and porosity
figure;
plotCellData(G, log10(rock.perm(:, 1)));
axis equal tight off
daspect([1 1 0.2])
view(85, 45); myZoom(1.2);
colorbar
title('Horizontal permeability (log10)')

figure;
plotCellData(G, log10(rock.perm(:, 3)));
axis equal tight off
daspect([1 1 0.2])
view(85, 45); myZoom(1.2);
colorbar
title('Vertical permeability (log10)')


figure;
plotCellData(G, rock.poro);
axis equal tight off
daspect([1 1 0.2])
view(85, 45); myZoom(1.2);
colorbar
title('Porosity')

%% Set up wells 
% We set up a number of somewhat arbitrary wells around the domain. We
% drain a complete pore volume using three producers, with four injectors
% at fixed bottom hole pressures that give pressure support.
totTime = 100*year;
N_step = 100;
dt = totTime/N_step;

pv = poreVolume(G, rock);

wells = [13, 88,  -1; ...
         18, 87,  -1; ...
         36, 90,  -1; ...
         10, 15,  -1; ...
         24, 32,  1; ...
         8,  45,  1; ...
         16, 55,  1];

W = [];
[inum, pnum] = deal(1);
for i = 1:size(wells, 1)
    % Set well
    W = verticalWell(W, G, rock, wells(i, 1), wells(i, 2), [],...
                     'comp_i', [1, 0], 'type', 'bhp');
    if wells(i, 3) == 1
        % Producer
        W(i).val = -sum(pv)/(totTime*sum(wells(:, 3) == 1));
        W(i).type = 'rate';

        W(i).name = ['P', num2str(pnum)];
        W(i).sign = -1;
        pnum = pnum + 1;
    else
        % Injector
        W(i).val = 500*barsa;
        W(i).sign = 1;
        W(i).name = ['I', num2str(inum)];
        inum = inum + 1;
    end
end

% Plot the grid, the wells and the perforated cells
close all
plotGrid(G, 'FaceColor', 'none', 'EdgeA', .2)
plotWell(G, W)
plotGrid(G, vertcat(W.cells), 'FaceColor', 'none', 'EdgeColor', 'b')
axis equal tight off
daspect([1 1 0.2])
view(80, 65); myZoom(1.5);

%% Simulate the base case
T = getFaceTransmissibility(G, rock);
fluid = initSimpleFluid('mu', [1, 5]*centi*poise, 'n', [2, 2], 'rho', [0, 0]);

state0 = initResSol(G, 0, [0, 1]);
gravity reset off

psolve = @(state) incompTPFA(state, G, T, fluid, 'Wells', W, 'use_trans', true);
solver = @(state) implicitTransport(state, G, dt, rock, fluid, 'wells', W);

state = psolve(state0);

states([1 N_step+1]) = state;
for i = 1:N_step
    fprintf('Step %d of %d: ', i, N_step);
    
    fprintf('Solving pressure... ');
    state = psolve(states(end));
    fprintf('Ok! ');
    fprintf('Solving transport... ');
    states(i+1) = solver(state);
    fprintf('Ok!\n');
end

%% Set up coarse grid
% We have two options. The first is to use a fully unstructured coarse grid
% based on Metis. The second is to simply use a combination of uniform
% partitioning and splitting of carse blocks over faults.
%
% We use METIS if available. If you have METIS installed, but are unable to
% use it in MRST, please see the documentation in callMetisMatrix.
global METISPATH
useMETIS = ~isempty(METISPATH);

cdims = ceil(G.cartDims./[15, 10, 10]);
if useMETIS
    p = partitionMETIS(G, T, 250);
    p = processPartition(G, p);
else
    padded = partitionUniformPadded(G, [cdims(1:2), 1]);
    uni = partitionUI(G, [1, 1, cdims(3)]);
    p = padded.*uni;

    G_fault = makeInternalBoundary(G, find(G.faces.tag > 0));
    p = processPartition(G_fault, p);
end
p = compressPartition(p);

% Merge smaller blocks
mrstModule add agglom
p0 = p;
fconn = ones(G.faces.num, 1);
fconn(G.faces.tag > 0) = 0;
p = mergeBlocksByConnections(G, p, fconn, 25);

p = processPartition(G, p);
p = compressPartition(p);

% Plot the coarse grid
CG = generateCoarseGrid(G, p);
figure; plotCellData(G, mod(p, 13), 'EdgeColor', 'none')
plotGrid(CG, 'facec', 'none', 'edgec', 'w', 'linewidth', 1)
axis equal tight off
daspect([1 1 0.2])
view(85, 45); myZoom(1.5);

%% Set up support regions required for MsRSB and move center points to wells
CG = coarsenGeometry(CG);
CG = addCoarseCenterPoints(CG);
CG = setCentersByWells(CG, W);
CG = storeInteractionRegion(CG, 'ensureConnected', true);
CG = setupMexInteractionMapping(CG);

%% Set up basis functions
% By default, we use the C-accelerated version 
useCompiledBasis = true;
A = getIncomp1PhMatrix(G, T);
% Update basis functions every now and then
updateBasis = true;

getBasis = @(A) getMultiscaleBasis(CG, A, 'type', 'MsRSB', 'useMex', useCompiledBasis);
basis0 = getBasis(A);

%% Solve multiscale without iterations
% We solve the base case where only the multiscale solver is used.
basis = basis0;
W_ms = W;
fn = getSmootherFunction('type', 'ilu');

psolve = @(state, basis) incompMultiscale(state, CG, T, fluid, basis, 'wells', W_ms, ...
    'getSmoother', fn, 'iterations', 0);

states_ms([1 N_step+1]) = psolve(state0, basis);
for i = 1:N_step
    state = states_ms(end);

    if updateBasis && mod(i, 10) == 0 && i > 1
        A = getIncomp1PhMatrix(G, T, state, fluid);
        basis = getBasis(A);
    end
    
    fprintf('Step %d of %d: ', i, N_step);
    fprintf('Solving pressure... ');
    state = psolve(state, basis);
    fprintf('Ok! ');
    fprintf('Solving transport... ');
    states_ms(i+1) = solver(state);
    fprintf('Ok!\n');
end

%% Solve using multiscale with additional iterations
% We can apply 5 multiscale-ILU(0) cycles at each step in order to improve
% the solution quality at a low cost. The pressure is still far from
% converged, but this will help the solver to get the impact of wells
% correctly.
psolve = @(state, basis) incompMultiscale(state, CG, T, fluid, basis, 'wells', W_ms, ...
    'getSmoother', fn, 'iterations', 5, 'useGMRES', true);
basis = basis0;

states_it([1 N_step+1]) = psolve(state0, basis);
for i = 1:N_step
    state = states_it(end);
    if updateBasis && mod(i, 10) == 0 && i > 1
        A = getIncomp1PhMatrix(G, T, state, fluid);
        basis = getBasis(A);
    end
    
    fprintf('Step %d of %d: ', i, N_step);
    fprintf('Solving pressure... ');
    state = psolve(state, basis);
    fprintf('Ok! ');
    fprintf('Solving transport... ');
    states_it(i+1) = solver(state);
    fprintf('Ok!\n');
end

%% Set up interactive plotting of reservoir states
names = {'Finescale', 'MsRSB', 'MsRSB (5 cycles)'};

close all; plotToolbar(G, states);
axis equal tight off
daspect([1 1 0.2])
view(85, 20);
plotWell(G, W);
title(names{1});
colorbar('horiz')

figure; plotToolbar(G, states_ms);
axis equal tight off
daspect([1 1 0.2])
view(85, 20);
plotWell(G, W);
title(names{2});
colorbar('horiz')

figure; plotToolbar(G, states_it);
axis equal tight off
daspect([1 1 0.2])
view(85, 20);
plotWell(G, W);
title(names{3});
colorbar('horiz')

%% Launch interactive plotting of wells
Time = cumsum([0; repmat(dt, N_step, 1)]);

ws_ref = convertIncompWellSols(W, states, fluid);
ws_ms = convertIncompWellSols(W, states_ms, fluid);
ws_it = convertIncompWellSols(W, states_it, fluid);

ws = {ws_ref, ws_ms, ws_it};
shortname = {'ref', 'ms', 'it'};
plotWellSols(ws, Time, 'datasetnames', names)

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
