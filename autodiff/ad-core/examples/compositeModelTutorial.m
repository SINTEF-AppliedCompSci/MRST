%% CompositeModel tutorial
% This tutorial shows how to set up and simulate a composite model, i.e., a
% model built up of multiple submodels. The framework can for instance be
% used to combine models with different physics, and for coupled simulation
% of multiple reservoirs.

% PLEASE NOTE: The composite model functionality is under active
% development, and will likely change significantly between releases.

%% Load modules
mrstModule add test-suite spe10
mrstModule add ad-blackoil ad-core ad-props
mrstModule add mrst-gui

mrstVerbose on

%% Set up model models
% We consider two layers from the SPE10 model 2 benchmark: Layer 1, from
% the Tarbert formation, and Layer 85, from the Upper Ness formation

% Use test-suite to set up two subsets of layer 1 and layer 85
layer01 = TestCase('spe10_wo', 'layers',  1, 'J', 1:110, 'make2D', true);
layer85 = TestCase('spe10_wo', 'layers', 85, 'J', 1:110, 'make2D', true);

% Layer 1: Pick well in lower left corner and set to inject at fixed BHP
layer01.schedule.control.W       = layer01.schedule.control.W(1);
layer01.schedule.control.W.type  = 'bhp';
layer01.schedule.control.W.val   = layer01.state0.pressure(1) + 1500*barsa;
layer01.schedule.control.W.compi = [1,0];
layer01.schedule.control.W.sign  = 1;

% Layer 85: Pick well in upper right corner and set to produce at fixed BHP
layer85.schedule.control.W = layer85.schedule.control.W(4);
layer85.schedule.control.W.type = 'bhp';
layer85.schedule.control.W.val  = layer85.state0.pressure(1) - 200*barsa;
layer85.schedule.control.W.sign = -1;

%% Make composite model
% Different models can be combined using CompositeModel. This takes a cell
% array of models as input. The submodels are stored in a property
% `submodels`, with field names that by default equals the class of the
% submodel. In this case, the submodels are of the same class, and we use
% the optional input argument `names` to give the two models unique names
model = CompositeModel({layer01.model, layer85.model}, ...
    'names', {'Layer01', 'Layer85'}); % Submodel names

% The schedule has the same structure as a regular schedule, except for
% that each control has one field per submodel. The field name must equal
% the corresponding name in the submodel
schedule = layer01.schedule;
schedule.control = struct( ...
    'Layer01', layer01.schedule.control, ...
    'Layer85', layer85.schedule.control  ...
);

% States follows a similar structure, with one field per submodel.
state0 = struct( ...
    'Layer01', layer01.state0, ...
    'Layer85', layer85.state0  ...
);

%% Define couplings
% In this case, we want to couple the models so that everything entering
% the upper-right cell in layer 1 is extracted and injected into the
% lower-right cell in layer 85. This is specified with a state function
% SourceCouplingTerm, which inherits from the base class CouplingTerm. All
% state functions of this class has a struct `couplings` that specifies
% which subset and which equation the sources/sinks sould be instered into.

% We can specify this explcitly as follows:
cell01 = layer01.model.G.cells.num  ; % Upper-right corner
cell85 = layer85.model.G.cartDims(1); % Lower-right corner
couplings = { ...
    struct( ... % Specify sink in upper-right corner of layer 1
        'model'    , 'Layer01'         , ... % Model name
        'equations', {{'water', 'oil'}}, ... % Equations
        'subset'   , cell01            , ... % Upper-right corner
        'sign'     , -1                  ... % Sink
    ), ...
    struct( ... Specify source in lower-right corner of layer 1
        'model'    , 'Layer85'         , ... % Model name
        'equations', {{'water', 'oil'}}, ... % Equations
        'subset'   , cell85            , ... % Lower-right corner
        'sign'     , 1                   ... % Source
    ), ...
};
% The coupling structs can then be provided as an optional input argument
% to our coupling term state function
sct = SourceCouplingTerm(model, ...
    'Layer01'  , ...
    'Layer85'  , ...
    'couplings', couplings ...
); %#ok
   
% Alternatively, when the coupled models share the same equations, we can
% let `SourceCouplingTerm` set up the couplings for us.
sct = SourceCouplingTerm(model, ...
    'Layer01'  , ...
    'Layer85'  , ...
    'cellsA'   , cell01, ...
    'cellsB'   , cell85  ...
);

% We then set the coupling to the model
model = model.setCouplingTerm(sct, 'Source');

% Store the coupling term output to state for visualization later
model.OutputStateFunctions = model.CouplingTerms.getNamesOfStateFunctions();

%% Simulate case
problem = packSimulationProblem(state0, model, schedule, 'compositeModelTutorial');
simulatePackedProblem(problem, 'restartStep', 1);

%% Load results
[wellSols, states, reports] = getPackedSimulatorOutput(problem);

%% Visualize
close all

% Make figure
df = get(0, 'DefaultFigurePosition');
figure('Position', [df(1:2), 750, 800]);

% Plot water saturation in layer 1
subplot(4,2,1:2:5);
ax01 = plotCellData(layer01.model.G, states{1}.Layer01.s(:,1), 'EdgeColor', 'none');
layer01.plotWells();
layer01.setAxisProperties(gca); caxis([0.2, 0.8]);
title('Layer 1')

% Plot water saturation in layer 85
subplot(4,2,2:2:6);
ax85 = plotCellData(layer85.model.G, states{1}.Layer85.s(:,1), 'EdgeColor', 'none');
layer85.plotWells();
layer85.setAxisProperties(gca); caxis([0.2, 0.8]); 
title('Layer 85')

colormap(pink);

% Plot mass flux betweent the two models
subplot(4,2,7:8);
q = cellfun(@(state) horzcat(state.CouplingTerms.Source{:}), states, 'UniformOutput', false);
q = abs(cell2mat(q));
time = cumsum(schedule.step.val)/year;
axq = plot(time, q.*nan, 'LineWidth', 2);
axis tight; box on, grid on, xlim([0,time(end)]); ylim([0, max(q(:))]);
legend({'Water mass rate', 'Oil mass rate'}, ...
    'Location', 'northwest');
xlabel('Time (years)'); ylabel('Rate (kg/s)');

for i = 1:numel(states)
    
    ax01.CData = states{i}.Layer01.s(:,1);
    ax85.CData = states{i}.Layer85.s(:,1);
    axq(1).YData(1:i) = q(1:i,1);
    axq(2).YData(1:i) = q(1:i,2);
    
    pause(0.1);
    
end 

%% Interacitve plot
layer01.plot(states)

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
