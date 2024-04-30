%% Calibrate Sleipner model 
% This example calibrates the Sleipner model (top surface elevations, rock
% property, CO2 entry rates, CO2 density) by minimizing the misfit between
% simulated and observed CO2 plume thicknesses.
%
% The calibration result obtained by running this script is only one
% possible solution that is found within a "family of optimals".
%
% For more details, see Nilsen et al. "Using sensitivities and
% vertical-equilibrium models for parameter estimation of CO2 injection
% models with application to Sleipner data", Energy Procedia, 2017.
% (Presented at the 13th International Conference on Greenhouse Gas Control
% and Technologies, GHGT-13, 14-18 November 2016, Lausanne, Switzerland)

mrstModule add co2lab-common co2lab-ve optimization
gravity('reset',[0 0 9.8])
gravity on

%% Set up model with a schedule of 12 years (1999-2010).
% An initial simulation will be run.
num_years = 12;
[ model, schedule, initState, states, wellSols ] = ...
    calibrateSleipnerSetup( 'num_years',num_years );

Gt = model.G;

schedule_base = schedule;
states_base   = states;
wellSols_base = wellSols;

%% Get 2001, 2004, 2006, 2010 plume heights data
assert(num_years == 12);
assert(numel(states) == 12);
clear newplumes
newplumes=cell(numel(states),1);
for ty=[3,6,8,12]
    newplumes{ty}.year = ty + 1998;
    newplumes{ty}.FH = getSleipnerPlumeHeights('year',newplumes{ty}.year);
    newplumes{ty}.h = newplumes{ty}.FH(Gt.cells.centroids(:,1), Gt.cells.centroids(:,2));
end

% Initial simulated heights vs. observed plume thicknesses
plotObsAndSim(Gt, newplumes, states_base);

%% Define parameters to optimize.
% The limits for grid elevation change (dzLims) refer to the
% dimensional unit of meters for which the initial grid can change. The
% limits for CO2 density, permeability, and porosity are relative 
% multipliers applied to the initial variable values. The limits for CO2
% entry rates refer to the (constant) multiplier applied to the
% initial 1999-2010 entry rates.

setup = struct('model', model, 'schedule', schedule, 'state0', initState);

params = [];

params = addParameter(params, setup, 'name', 'dz', ...
                      'belongsTo', 'model', ...
                      'location', {'G', 'dz'}, ...
                      'boxLims', [-5, 5], 'uniformLimits', true);

params = addParameter(params, setup, 'name', 'perm', ...
                      'type', 'multiplier', ...
                      'belongsTo', 'model', ...
                      'location', {'rock', 'perm'}, ...
                      'setfun', @setIsotropicPermeabilityFun, ...
                      'lumping', ones(Gt.cells.num, 1), ...
                      'relativeLimits', [0.5 12], 'scaling', 'linear');

params = addParameter(params, setup, 'name', 'porevolume', ...
                      'type', 'multiplier', ...
                      'lumping', ones(Gt.cells.num, 1), ...
                      'relativeLimits', [0.5 1.4]);

params = addParameter(params, setup, 'name', 'rhoGmult', ...
                      'type', 'multiplier', ...
                      'belongsTo', 'model', ...
                      'location', {'fluid', 'rhoGmult'}, ...
                      'setfun', @setrhoGmultfun, ...
                      'getfun', @getrhoGmultfun, ...
                      'relativeLimits' , [0.4, 2], 'scaling', 'linear');

for i = 1:numel(setup.schedule.control)
    params = addParameter(params, setup, 'name', 'rate', ...
                          'type', 'value', ...
                          'controlSteps', i, ...
                          'relativeLimits', [0.1 1.5]); 
end


%% Define objective function and compute initial value

% Determine scaling by comparing with base result
v_base = matchToPlumeData(model, states_base, newplumes);
scaling = max(sum(abs(cell2mat(v_base))), 1);

% Initial parameter values
p0 = getScaledParameterVector(setup, params);

% Define an objective function p -> obj(p) that can be used in the 
% nonlinear optimization routine
objh = @(p) evaluatePlumeMatch(p, @matchToPlumeData, setup, params, newplumes, scaling);
                               

% Run optimization (parameter calibration)
[v, p_opt, history] = unitBoxBFGS(p0, objh);

% Recompute optimal solution
[v_opt, der_opt, wellSols_opt, states_opt, setup_opt] = objh(p_opt);

% Plot final simulated heights vs. observed plume thicknesses
states_opt = addHeightData(states_opt, Gt, model.fluid); % for easier plotting
plotObsAndSim(Gt, newplumes, states_opt)


%% Some final output
cvols = setup.model.G.cells.volumes .* setup.model.G.cells.H;
q_mean = mean(arrayfun(@(x)x.W.val, setup.schedule.control));
q_mean_opt = mean(arrayfun(@(x)x.W.val, setup_opt.schedule.control));

fprintf('\nVariable             Initial    Calibrated    Unit \n')
fprintf('----------------------------------------------------\n')
fprintf(' avg. porosity      |  %1.3f   |  %1.3f      |      \n', ...
        mean(setup.model.operators.pv ./ cvols), mean(setup_opt.model.operators.pv ./ cvols));
fprintf(' avg. permeability  |  %2.3f   |  %2.3f     | darcy\n', ...
        convertTo(mean(setup.model.rock.perm), darcy), convertTo(mean(setup_opt.model.rock.perm),darcy));
fprintf(' CO2 density        |  %4.1f   |  %4.1f      | kg/m3\n', ...
        setup.model.fluid.rhoGS, setup_opt.model.fluid.rhoGS);
fprintf(' avg. CO2 entry rate| %0.4f   | %0.4f      | m3/s \n', q_mean, q_mean_opt);

% Initial versus calibrated grid:
figure
subplot(1,3,1)
title({'Initial';'top-surface elevations'})
plotCellData(setup.model.G, setup.model.G.cells.z, 'edgealpha',0.1); colorbar; axis equal tight off
subplot(1,3,2)
title({'Calibrated';'top-surface elevations'})
plotCellData(setup_opt.model.G, setup_opt.model.G.cells.z + setup_opt.model.G.dz, 'edgealpha',0.1); colorbar; axis equal tight off
subplot(1,3,3)
title({'Difference'; 'top-surface elevations'});
plotCellData(setup_opt.model.G, setup_opt.model.G.dz, 'edgealpha',0.1); colorbar; axis equal tight off

%% Copyright Notice
%
% <html>
% <p><font size="-1">
% Copyright 2009-2024 SINTEF Digital, Mathematics & Cybernetics.
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
