% sensitivitiesModel2D - analyse sensitivity capabilities 

mrstModule add ad-core ad-blackoil ad-props optimization spe10 mrst-gui

% setup model -> grid, rock, schedule, fluid etc
setupModel2D

% Create model-object of class TwoPhaseOilWaterModel
model  = TwoPhaseOilWaterModel(G, rock, fluid);
% Set initial state and run simulation:
state0 = initResSol(G, 200*barsa, [0, 1]);

% Set up a reference ("truth") model with different pv and trans:
model_ref = model;
model_ref.operators.pv = model.operators.pv + 50*ones(G.cells.num,1);
model_ref.operators.T  = model.operators.T.*(1-.3*rand(numel(model.operators.T),1));

% run ref model
[ws_ref, states_ref, r_ref] = simulateScheduleAD(state0, model_ref, schedule);
% run model
[ws, states, r] = simulateScheduleAD(state0, model, schedule);

% plot well solutions for the two models
plotWellSols({ws_ref, ws}, {r_ref.ReservoirTime, r.ReservoirTime}, ...
            'datasetnames', {'reference', 'current'})


        %% setup misfit-function and run adjoint to get parameter sensitivities
% choose parameters (trans and pv currently only opotions)
params      = {'transmissibility', 'porevolume'}; 
% treat actual property value as parameter (se below for multiplyers)
paramTypes  = {'value', 'value'};
% setup weights for matching function, empty weight uses default (will 
% produce weighting of ~O(1)). 
weighting =  {'WaterRateWeight',     [], ...
              'OilRateWeight',       [] , ...
              'BHPWeight',           []};
   
% compute misfit function value (first each summand corresonding to each time-step)
misfitVals = matchObservedOW(G, ws, schedule, ws_ref, weighting{:});
% sum values to obtiain scalar objective 
misfitVal = sum(vertcat(misfitVals{:}));
fprintf('Current misfit value: %6.4f\n', misfitVal)

% setup (per time step) mismatch function handle for passing on to adjoint sim
objh = @(tstep)matchObservedOW(G, ws, schedule, ws_ref, 'computePartials', true, 'tstep', tstep, weighting{:});

% run adjoint to compute sensitivities of misfit wrt params, sensitivities
% are output in same order as parameters are listed
sens = computeSensitivitiesAdjointAD(state0, states, model, schedule, objh, ...
                                     'Parameters'    , params, ...
                                     'ParameterTypes', paramTypes);

%% plot sensitivities on grid:
figure,
subplot(1,3,1), plotCellData(G, log(rock.perm(:,1)), 'EdgeColor', 'none'), title('log permeability')
plotWellData(G, W);colorbar
subplot(1,3,2), plotCellData(G, sens{2}, 'EdgeColor', 'none'), colorbar,title('Pore volume sensitivity');
subplot(1,3,3), plotCellData(G,  cellAverage(G, sens{1}), 'EdgeColor', 'none'), colorbar,title('Average transmissiblity sensitivity');


%% do the same exercise as above with multiplyers rather than values 
% note that the gradient wrt (unit) multiplyer is equal to the gradient wrt value devided by value)
paramTypes  = {'multiplyer', 'multiplyer'};
sens = computeSensitivitiesAdjointAD(state0, states, model, schedule, objh, ...
                                     'Parameters'    , params, ...
                                     'ParameterTypes', paramTypes);

%% Same plotting as above but now sensitivities wrt multiplyers:
figure,
subplot(1,3,1), plotCellData(G, log(rock.perm(:,1)), 'EdgeColor', 'none'), title('log permeability')
plotWellData(G, W);colorbar
subplot(1,3,2), plotCellData(G, sens{2}, 'EdgeColor', 'none'), colorbar,title('Pore volume multiplyer sensitivity');
subplot(1,3,3), plotCellData(G,  cellAverage(G, sens{1}), 'EdgeColor', 'none'), colorbar,title('Average transmissiblity multiplyer sensitivity');








