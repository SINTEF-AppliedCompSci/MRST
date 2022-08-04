%% 3D, two-phase simulation of the Johansen formation
% This example shows the creation of a 3D grid for the Johansen formation using
% the inbuilt mrst utility function makeJohansenVEgrid. This utility function
% creates both the 3D and the VE grid for Johansen but in this case we only use
% the 3D grid. We then run a two-phase simulation with 100 years of CO2 injection
% into a brine filled reservoir followed by 1000 years post injection modelling.
% This example is similar to basic_example_3D.m but with a more realistic grid.

%% Load modules

mrstModule add ad-core;
mrstModule add ad-props;
mrstModule add ad-blackoil;

%% Grid and rock

% Load Johansen model
[G, rock, bcIx] = makeJohansenVEgrid();

% Plot grid and rock properties
figure; plotCellData(G, rock.poro); view(-35,15); colorbar; % porosity
set(gcf, 'position', [531   337   923   356]); axis tight;
figure; plotCellData(G, rock.perm(:,1)/darcy); view(-55, 60); colorbar; % lateral permeability
set(gcf, 'position', [152   419   1846   700], 'color', 'white'); axis tight;
set(gca, 'fontsize', 24);
figure; plotCellData(G, rock.perm(:,3)/darcy); view(-55, 60); colorbar; % vertical permeability
set(gcf, 'position', [152   419   1846   700], 'color', 'white'); axis tight;
set(gca, 'fontsize', 24);

%% Initial state
gravity on; % tell MRST to turn on gravity
g = gravity; % get the gravity vector
rhow = 1000; % density of brine corresponding to 94 degrees C and 300 bar
initState.pressure = rhow * g(3) * G.cells.centroids(:,3);
initState.s = repmat([1, 0], G.cells.num, 1);
initState.sGmax = initState.s(:,2);

%% Fluid model
co2     = CO2props(); % load sampled tables of co2 fluid properties
p_ref   = 30 * mega * Pascal; % choose reference pressure
t_ref   = 94 + 273.15; % choose reference temperature, in Kelvin
rhoc    = co2.rho(p_ref, t_ref); % co2 density at ref. press/temp
cf_co2  = co2.rhoDP(p_ref, t_ref) / rhoc; % co2 compressibility
cf_wat  = 0; % brine compressibility (zero)
cf_rock = 4.35e-5 / barsa; % rock compressibility
muw     = 8e-4 * Pascal * second; % brine viscosity
muco2   = co2.mu(p_ref, t_ref) * Pascal * second; % co2 viscosity

mrstModule add ad-props; % The module where initSimpleADIFluid is found

% Use function 'initSimpleADIFluid' to make a simple fluid object
fluid = initSimpleADIFluid('phases', 'WG'           , ...
                           'mu'  , [muw, muco2]     , ...
                           'rho' , [rhow, rhoc]     , ...
                           'pRef', p_ref            , ...
                           'c'   , [cf_wat, cf_co2] , ...
                           'cR'  , cf_rock          , ...
                           'n'   , [2 2]);
% plot relperm curves
sw = linspace(0, 1, 200);
figure; hold on;
plot(sw, fluid.krW(sw),   'b', 'linewidth', 1.5);
plot(sw, fluid.krG(1-sw), 'r', 'linewidth', 1.5);
xlabel('brine saturation'); ylabel('relative permeability')
set(gca, 'fontsize', 14);

% Change relperm curves
srw = 0.27;
src = 0.20;
fluid.krW = @(s) fluid.krW(max((s-srw)./(1-srw), 0));
fluid.krG = @(s) fluid.krG(max((s-src)./(1-src), 0));

% Plot relperm curves
figure; hold on;
sw = linspace(srw, 1, 200);
plot(sw, fluid.krW(sw),   'b', 'linewidth', 1.5);
plot(sw, fluid.krG(1-sw), 'r', 'linewidth', 1.5);
line([srw, srw], [0 1], 'color', 'black', 'linestyle', ':', 'linewidth', 1);
xlabel('brine saturation'); ylabel('relative permeability')
set(gca, 'fontsize', 14, 'xlim', [0 1]);


% Add capillary pressure curve
pe = 5 * kilo * Pascal;
pcWG = @(sw) pe * sw.^(-1/2);
fluid.pcWG = @(sg) pcWG(max((1-sg-srw)./(1-srw), 1e-5)); %@@


%% Wells

% Well cell indices in 'global' grid: 48, 48, 6:10
wc_global = false(G.cartDims);
wc_global(48, 48, 6:10) = true;
wc = find(wc_global(G.cells.indexMap));

% plot well cells
plotGrid(G, 'facecolor', 'none', 'edgealpha', 0.1);
plotGrid(G, wc, 'facecolor', 'red');


% Calculate the injection rate
inj_rate = 3.5 * mega * 1e3 / year / fluid.rhoGS;

% Start with empty set of wells
W = [];

% Add a well to the set
W = addWell([], G, rock, wc, ...
            'type', 'rate', ...  % inject at constant rate
            'val', inj_rate, ... % volumetric injection rate
            'comp_i', [0 1]);    % inject CO2, not water

%% Boundary conditions

% Start with an empty set of boundary faces
bc = [];

% Computing hydrostatic pressure for boundary cells
p_bc = G.faces.centroids(bcIx, 3) * rhow * g(3);

% Add hydrostatic pressure conditions to open boundary faces
bc = addBC(bc, bcIx, 'pressure', p_bc, 'sat', [1, 0]);

%% Schedule

% Setting up two copies of the well and boundary specifications. 
% Modifying the well in the second copy to have a zero flow rate.
schedule.control    = struct('W', W, 'bc', bc);
schedule.control(2) = struct('W', W, 'bc', bc);
schedule.control(2).W.val = 0;

% Specifying length of simulation timesteps
schedule.step.val = [repmat(year, 100, 1); ...
                     repmat(10*year, 100, 1)];

% Specifying which control to use for each timestep.
% The first 100 timesteps will use control 1, the last 100
% timesteps will use control 2.
schedule.step.control = [ones(100, 1); ...
                         ones(100, 1) * 2];

%% Model

model = TwoPhaseWaterGasModel(G, rock, fluid, 0, 0);

%% Simulate
[wellSol, states] = simulateScheduleAD(initState, model, schedule);

%% Plot results
figure; plotCellData(G, states{100}.s(:,2)); view(-63, 68); colorbar; 
set(gcf, 'position', [531   337   923   356]); axis tight;
figure; plotCellData(G, states{200}.s(:,2)); view(-63, 68); colorbar; 
set(gcf, 'position', [531   337   923   356]); axis tight;

% vertical plot
[i j k] = ind2sub(G.cartDims, G.cells.indexMap);
figure; plotCellData(extractSubgrid(G, j==48 & i>18 & i < 60), states{100}.s(j==48 & i>18 & i < 60,2)); 
view(0,0); axis tight;

figure; plotCellData(extractSubgrid(G, j==48 & i>18 & i < 60), states{200}.s(j==48 & i>18 & i < 60,2)); 
view(0,0); axis tight;

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

