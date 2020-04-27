%% Basic VE example
% This example shows how to setup a basic VE simulation in mrst. We use the
% Johansen formation and create a VE grid from a 3D mrst grid. We demonstrate the
% use of an EOS and a model with a capillary fringe. We also analyse the
% evolution of trapping mechanisms throughout the simulation.
%
% This example is described in more detail in section 3.1 of
% "Simplified models for numerical simulation of geological CO2 storage"
% (O. Andersen, 2017)
%
% For comparison example3D.m shows the fully 3D version of this model.

%% Load modules

mrstModule add ad-core;
mrstModule add ad-props;

%% Setup model parameters

gravity on;
g       = gravity;
rhow    = 1000;
co2     = CO2props();
p_ref   = 30 *mega*Pascal;
t_ref   = 94+273.15;
co2_rho = co2.rho(p_ref, t_ref);
co2_c   = co2.rhoDP(p_ref, t_ref) / co2_rho;
wat_c   = 0; 
c_rock  = 4.35e-5 / barsa;
srw     = 0.27; % residual water
src     = 0.20; % residual CO2
pe      = 5 * kilo * Pascal; % capillary entry pressure
muw     = 8e-4 * Pascal * second; % brine viscosity
muco2   = co2.mu(p_ref, t_ref) * Pascal * second; % co2 viscosity

% Load Johansen model
[G, rock, bcIx, ~, ~, bcIxVE] = makeJohansenVEgrid();

% Specify well information
wc_global = false(G.cartDims);
wc_global(48, 48, 6:10) = true;
wc = find(wc_global(G.cells.indexMap));
%wc = [3715, 10210, 16022, 21396, 26770]';
inj_rate = 3.5 * mega * 1e3 / year / co2_rho;
W = addWell([], G, rock, wc, ...
            'type', 'rate', ...  % inject at constant rate
            'val', inj_rate, ... % volumetric injection rate
            'comp_i', [0 1]);    % inject CO2, not water

%% Grid and rock
% Make top surface grid used for VE simulation
[Gt, G, transMult] = topSurfaceGrid(G);

% Shift G 100 meters down, and plot both grids for comparison
GG = G;
GG.nodes.coords(:,3) = GG.nodes.coords(:,3) + 100;
figure; 
plotGrid(GG, 'facecolor', 'green'); 
plotGrid(Gt, 'facecolor', 'red');
view(-65,33);
set(gcf, 'position', [531 337 923 356]); axis tight;

% Computing vertically averaged rock object
rock2D = averageRock(rock, Gt);

%% Initial state
% Gt.cells.z gives the caprock depth of each cell in the 2D grid.  
initState.pressure = rhow * g(3) * Gt.cells.z;
initState.s = repmat([1, 0], Gt.cells.num, 1);
initState.sGmax = initState.s(:,2);

%% Fluid model

invPc3D = @(pc) (1-srw) .* (pe./max(pc, pe)).^2 + srw;
kr3D = @(s) max((s-src)./(1-src), 0).^2; % uses CO2 saturation

% Use full EOS rather than compressibilities and include capillary fringe model
% based on upscaled sampled capillary pressure.
% See makeVEFluid.m for a description of various fluid models which are
% available.
% A description of the different models can be found in the paper
% "Fully-Implicit Simulation of Vertical-Equilibrium Models with Hysteresis and
% Capillary Fringe" (Nilsen et al., Computational Geosciences 20, 2016).

fluid = makeVEFluid(Gt, rock, 'P-scaled table'             , ...
                    'co2_mu_ref'  , muco2, ...%6e-5 * Pascal * second , ...
                    'wat_mu_ref'  , muw, ...%8e-4 * Pascal * second , ...
                    'co2_rho_ref' , co2_rho                , ...
                    'wat_rho_ref' , rhow                   , ...
                    'co2_rho_pvt' , [co2_c, p_ref]         , ...
                    'wat_rho_pvt' , [wat_c, p_ref]         , ...
                    'residual'    , [srw, src]             , ...
                    'pvMult_p_ref', p_ref                  , ...
                    'pvMult_fac'  , c_rock                 , ...
                    'invPc3D'     , invPc3D                , ...
                    'kr3D'        , kr3D                   , ...
                    'transMult'   , transMult);

%% Set up schedule
W2D  = convertwellsVE(W, G, Gt, rock2D);

% hydrostatic pressure conditions for open boundary faces
p_bc = Gt.faces.z(bcIxVE) * rhow * g(3);
bc2D = addBC([], bcIxVE, 'pressure', p_bc); 
bc2D.sat = repmat([1 0], numel(bcIxVE), 1);

% Setting up two copies of the well and boundary specifications. 
% Modifying the well in the second copy to have a zero flow rate.
schedule.control    = struct('W', W2D, 'bc', bc2D);
schedule.control(2) = struct('W', W2D, 'bc', bc2D);
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
model = CO2VEBlackOilTypeModel(Gt, rock2D, fluid);

%% Simulate
[wellSol, states] = simulateScheduleAD(initState, model, schedule);

%% Plot results
% Plotting CO2 saturation for timestep 100 (100 years after start)
clf;
[h, h_max] = upscaledSat2height(states{100}.s(:,2), states{100}.sGmax, Gt, ...
                                'pcWG', fluid.pcWG, ...
                                'rhoW', fluid.rhoW, ...
                                'rhoG', fluid.rhoG, ...
                                'p', states{100}.pressure);
plotCellData(Gt.parent, height2Sat(struct('h', h, 'h_max', h_max), Gt, fluid));
colorbar; view(-63, 68); set(gcf, 'position', [531   337   923   356]); axis tight; 


% Plotting CO2 staturation for timestep 200 (1100 years after start)
[h, h_max] = upscaledSat2height(states{end}.s(:,2), states{end}.sGmax, Gt, ...
                                'pcWG', fluid.pcWG, ...
                                'rhoW', fluid.rhoW, ...
                                'rhoG', fluid.rhoG, ...
                                'p', states{end}.pressure);
plotCellData(Gt.parent, height2Sat(struct('h', h, 'h_max', h_max), Gt, fluid));
colorbar; view(-63, 68); set(gcf, 'position', [531   337   923   356]); axis tight; 

%% Plot trapping inventory
ta = trapAnalysis(Gt, false);
reports = makeReports(Gt, {initState states{:}}, model.rock, model.fluid, ...
                      schedule, [srw, src], ta, []);

h1 = figure; plot(1); ax = get(h1, 'currentaxes');
plotTrappingDistribution(ax, reports, 'legend_location', 'northwest');
set(gcf, 'position', [0 0 1100, 740])
set(gca, 'fontsize', 20);

%%
% <html>
% <p><font size="-1">
% Copyright 2009-2020 SINTEF Digital, Mathematics & Cybernetics.
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

