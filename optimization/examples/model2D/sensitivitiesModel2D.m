%% sensitivitiesModel2D - analyse sensitivity capabilities 

mrstModule add ad-core ad-blackoil ad-props optimization spe10 mrst-gui deckformat

% Setup model -> grid, rock, schedule, fluid etc
setupModel2D
%% Reset fluid to include scaling:
fluid = initSimpleADIFluid('mu',    [.3, 5, 0]*centi*poise, ...
                           'rho',   [1000, 700, 0]*kilogram/meter^3, ...
                           'n',     [2, 2, 0]);
                             
% saturation scaling
fluid.krPts  = struct('w', [0 0 1 1], 'ow', [0 0 1 1]);

scaling = {'SWL', .1, 'SWCR', .2, 'SWU', .9, 'SOWCR', .1, 'KRW', .9, 'KRO', .8};
% Create model-object of class TwoPhaseOilWaterModel
model_ref = TwoPhaseOilWaterModel(G, rock, fluid);                       
model_ref = imposeRelpermScaling(model_ref, scaling{:});


% Set initial state and run simulation:
state0 = initResSol(G, 200*barsa, [0, 1]); 

% Set up a perturbed model with different pv and perm:
rock1 = rock;
rock1.perm = rock.perm*1.1;
model =  GenericBlackOilModel(G, rock1, fluid, 'gas', false);  % TODO add a warning to only work with generic black oil
model = imposeRelpermScaling(model, scaling{:});
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
misfitVals = matchObservedOW(model, states, schedule, states_ref, weighting{:});
% sum values to obtiain scalar objective 
misfitVal = sum(vertcat(misfitVals{:}));
fprintf('Current misfit value: %6.4e\n', misfitVal)

% setup (per time step) mismatch function handle for passing on to adjoint sim
objh = @(tstep,model, state) matchObservedOW(model, states, schedule, states_ref,...
    'computePartials', true, 'tstep', tstep, weighting{:},'state',state);

%% run adjoint to compute sensitivities of misfit wrt params
% choose parameters, get multiplier sensitivities except for endpoints

model.toleranceCNV = 1e-6;
SimulatorSetup = struct('model', model, 'schedule', schedule, 'state0', state0);
parameters = [];
parameters = addParameter(parameters, SimulatorSetup, 'name', 'porevolume','type','multiplier');
parameters = addParameter(parameters, SimulatorSetup, 'name', 'permx','type','multiplier');
parameters = addParameter(parameters, SimulatorSetup, 'name', 'permy','type','multiplier');
parameters = addParameter(parameters, SimulatorSetup, 'name', 'permz','type','multiplier');
nms_relperm = {'swcr',  'sowcr', 'kro', 'krw'};
for k = 1:numel(nms_relperm)
   parameters = addParameter(parameters, SimulatorSetup, 'name', nms_relperm{k}); %, 'lumping',ones(n_cells,1)
end

raw_sens = computeSensitivitiesAdjointAD(SimulatorSetup, states,parameters, objh);

% do scaling of gradient for Multipliers
for kp = 1:numel(parameters)
    nm = parameters{kp}.name;
    if strcmp(parameters{kp}.type, 'multiplier')
        sens.(nm) = (raw_sens.(nm)).*parameters{kp}.getParameterValue(SimulatorSetup);
    else
        sens.(nm) = raw_sens.(nm);
    end
end
    
    
    
%% Plot sensitivities on grid:
figure,
subplot(2,2,1), plotCellData(G, log(rock.perm(:,1)), 'EdgeColor', 'none'), title('log permeability')
plotWellData(G, W);colorbar
subplot(2,2,2), plotCellData(G, sens.porevolume, 'EdgeColor', 'none'), colorbar,title('PV multiplier sensitivity');
subplot(2,2,3), plotCellData(G, sens.permx, 'EdgeColor', 'none'), colorbar,title('PermX multiplier sensitivity');
subplot(2,2,4), plotCellData(G, sens.permy, 'EdgeColor', 'none'), colorbar,title('PermY multiplier sensitivity');
%% Rel-perm end-point sensitivities
figure,
nms = {'swcr', 'sowcr', 'krw-max', 'kro-max'};
for k = 1:4
    subplot(2,2,k), plotCellData(G, sens.(parameters{k+2}.name), 'EdgeColor', 'none'), colorbar,title(nms{k});
end



%% Run new adjoint to obtain transmissibility and well connection transmissibilitiy sensitivities 

SimulatorSetup = struct('model', model, 'schedule', schedule, 'state0', state0);
parameters = [];
parameters = addParameter(parameters, SimulatorSetup, 'name', 'transmissibility','type','multiplier');
parameters = addParameter(parameters, SimulatorSetup, 'name', 'conntrans','type','multiplier');

raw_sens = computeSensitivitiesAdjointAD(SimulatorSetup, states,parameters, objh);
% do scaling of gradient for Multipliers
for kp = 1:numel(parameters)
    nm = parameters{kp}.name;
    if strcmp(parameters{kp}.type, 'multiplier')
        sens.(nm) = (raw_sens.(nm)).*parameters{kp}.getParameterValue(SimulatorSetup);
    else
        sens.(nm) = raw_sens.(nm);
    end
end

                                                                                             

%%
figure,
subplot(1,2,1), plotCellData(G,  cellAverage(G, sens.transmissibility), 'EdgeColor', 'none'), colorbar,title('Average trans multiplier sensitivity');
subplot(1,2,2), plot(sens.conntrans, '--o', 'MarkerSize', 14); title('Well connection trans multiplier sensitivity')
a = gca; a.XTick = 1:4;  a.XTickLabel = {W.name}; a.XLim; a.XLim = [.5 4.5];


%%
% <html>
% <p><font size="-1">
% Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.
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
