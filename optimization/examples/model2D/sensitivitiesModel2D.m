%% sensitivitiesModel2D - analyse sensitivity capabilities 

mrstModule add ad-core ad-blackoil ad-props optimization spe10 mrst-gui

% Setup model -> grid, rock, schedule, fluid etc
setupModel2D
%% Reset fluid to include scaling:
% $s_w -> \frac{s_w-swcr}{swu-swcr}$
% $s_o -> \frac{s_o-sowcr}{1-swl-sowcr}$
fluid = initSimpleScaledADIFluid('mu',    [.3, 5, 0]*centi*poise, ...
                                 'rho',   [1000, 700, 0]*kilogram/meter^3, ...
                                 'n',     [2, 2, 0], ...
                                 'swl',   0.10*ones(G.cells.num,1), ...
                                 'swcr',  0.15*ones(G.cells.num,1), ...
                                 'sowcr', 0.12*ones(G.cells.num,1), ...
                                 'swu',   0.90*ones(G.cells.num,1));
                                 
                       
% Create model-object of class TwoPhaseOilWaterModel
model_ref  = TwoPhaseOilWaterModel(G, rock, fluid, 'OutputStateFunctions', {});                       

% Set initial state and run simulation:
state0 = initResSol(G, 200*barsa, [.15, .85]); 

% Set up a perturbed model with different pv and perm:
rock1 = rock;
rock1.perm = rock.perm*1.1;
model = TwoPhaseOilWaterModel(G, rock1, fluid, 'OutputStateFunctions', {});       
model.operators.pv = model_ref.operators.pv.*0.8;

% run ref model
[ws_ref, states_ref, r_ref] = simulateScheduleAD(state0, model_ref, schedule);
% run model
[ws, states, r] = simulateScheduleAD(state0, model, schedule);

% plot well solutions for the two models
plotWellSols({ws_ref, ws}, {r_ref.ReservoirTime, r.ReservoirTime}, ...
            'datasetnames', {'reference', 'perturbed'})

%% setup misfit-function and run adjoint to get parameter sensitivities
% setup weights for matching function, empty weight uses default (will 
% produce function value of ~O(1) for 100% misfit). Only match rates in this example: 
weighting =  {'WaterRateWeight',     [], ...
              'OilRateWeight',       [] , ...
              'BHPWeight',           0};
   
% compute misfit function value (first each summand corresonding to each time-step)
misfitVals = matchObservedOW(G, ws, schedule, ws_ref, weighting{:});
% sum values to obtiain scalar objective 
misfitVal = sum(vertcat(misfitVals{:}));
fprintf('Current misfit value: %6.4e\n', misfitVal)

% setup (per time step) mismatch function handle for passing on to adjoint sim
objh = @(tstep)matchObservedOW(G, ws, schedule, ws_ref, 'computePartials', true, 'tstep', tstep, weighting{:});

% run adjoint to compute sensitivities of misfit wrt params
% choose parameters, get multiplier sensitivities except for endpoints
params      = {'porevolume', 'permeability', 'swl',   'swcr',  'sowcr', 'swu'}; 
paramTypes  = {'multiplier', 'multiplier',   'value', 'value', 'value', 'value'};   

sens = computeSensitivitiesAdjointAD(state0, states, model, schedule, objh, ...
                                     'Parameters'    , params, ...
                                     'ParameterTypes', paramTypes);

%% Plot sensitivities on grid:
figure,
subplot(2,2,1), plotCellData(G, log(rock.perm(:,1)), 'EdgeColor', 'none'), title('log permeability')
plotWellData(G, W);colorbar
subplot(2,2,2), plotCellData(G, sens.porevolume, 'EdgeColor', 'none'), colorbar,title('PV multiplier sensitivity');
subplot(2,2,3), plotCellData(G, sens.permx, 'EdgeColor', 'none'), colorbar,title('PermX multiplier sensitivity');
subplot(2,2,4), plotCellData(G, sens.permy, 'EdgeColor', 'none'), colorbar,title('PermY multiplier sensitivity');
%% Rel-perm end-point sensitivities
figure,
nms = {'Lower S_w', 'Critical S_w', 'Critical S_o', 'Upper S_w'};
for k = 1:4
    subplot(2,2,k), plotCellData(G, sens.(params{k+2}), 'EdgeColor', 'none'), colorbar,title(nms{k});
end


%% Run new adjoint to obtain transmissibility and well connection transmissibilitiy sensitivities 
params      = {'transmissibility', 'conntrans'};
paramTypes  = {'multiplier', 'multiplier'};
sens = computeSensitivitiesAdjointAD(state0, states, model, schedule, objh, ...
                                     'Parameters'    , params, ...
                                     'ParameterTypes', paramTypes);

figure,
subplot(1,2,1), plotCellData(G,  cellAverage(G, sens.transmissibility), 'EdgeColor', 'none'), colorbar,title('Average trans multiplier sensitivity');
subplot(1,2,2), plot(sens.conntrans, '--o', 'MarkerSize', 14); title('Well connection trans multiplier sensitivity')
a = gca; a.XTick = 1:4;  a.XTickLabel = {W.name}; a.XLim; a.XLim = [.5 4.5];


model.fluid = initSimpleScaledADIFluid(model.fluid, 'swl', zeros(G.cells.num,1), ...
                                                    'swcr', zeros(G.cells.num,1), ...
                                                    'sowcr', zeros(G.cells.num,1), ...
                                                    'swu', ones(G.cells.num,1));
sens = computeSensitivitiesAdjointAD(state0, states, model, schedule, objh, ...
                                     'Parameters'    , {'swl', 'swcr', 'sowcr', 'swu'});

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
