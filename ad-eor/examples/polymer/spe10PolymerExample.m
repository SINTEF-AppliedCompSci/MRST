%% Black-Oil Polymer System for a Layer of the SPE10 Model
%
mrstModule add ad-core ad-blackoil ad-eor ad-props deckformat mrst-gui spe10

%% Use setupSPE10_AD to Fetch the SPE10 model
% We pick up only one layer 
%
layers = 35;
[~, model, ~] = setupSPE10_AD('layers', layers);
% We recover the grid and rock properties from the model
G = model.G;
rock = model.rock;


%% Fluid Properties
%
% We use the same blackoil properties as in the 2D case.
%
fname = {'BOPOLYMER.DATA', ...
         'POLY.inc', ...
         'BOPOLYMER_NOSHEAR.DATA', ...
         'POLY_NOSHEAR.inc', ...
        };
files = fullfile(getDatasetPath('BlackoilPolymer2D', 'download', true), fname);

% check to make sure the files are complete
e = cellfun(@(pth) exist(pth, 'file') == 2, files);

if ~all(e)
    pl = ''; if sum(e) ~= 1, pl = 's'; end
    msg = sprintf('Missing data file%s\n', pl);
    msg = [msg, sprintf('  * %s\n', fname{~e})];
    error('Dataset:Incomplete', msg);
end


% Parsing the data file with shear-thinning effect.
deck = readEclipseDeck(files{1});
% The deck is using metric system, MRST uses SI unit internally
deck = convertDeckUnits(deck);

% fluid properties, such as densities, viscosities, relative
% permeability, etc.
fluid = initDeckADIFluid(deck);

%% Constructing the physical model used for this simulation
% Here, we use three phase blackoil-polymer model
model = ThreePhaseBlackOilPolymerModel(G, rock, fluid);
model.disgas = true;
model.vapoil = true;


%% Define the Wells
% The wells are set to operate with a high rate and no pressure limit.
% Hence, the pressure will rise to far above what can be used in real
% operational settings. However, the main point of the example is to force
% flow through large parts of the heterogeneous reservoir and observe the
% evolution of the displacement fronts.
injeIJ = [59  17];        % Location of injection well
prodIJ = [ 2 194];        % Location of production well
rate   = 2*meter^3/day;   % Injection rate
bhp    = 200*barsa;       % Pressure at production well
nz     = G.cartDims(3);

W = [];
% Set up injection well (rate control)
% The polymer injection concentration is set later, see below
W = verticalWell(W, G, rock, injeIJ(1), injeIJ(2), 1:nz, ...
                 'Type'   , 'rate', ...
                 'Val'    , rate, ...
                 'Radius' , 0.1, ...
                 'Comp_i' , [1, 0, 0], ...
                 'cp'     , 0, ...
                 'name'   , 'INJE', ...
                 'Sign'   , 1);
% Set up production well (pressure control)
bhpProd = 100*barsa;
W = verticalWell(W, G, rock, prodIJ(1), prodIJ(2), 1:nz, ...
                 'Type'   , 'bhp', ...
                 'Val'    , bhpProd, ...
                 'Radius' , 0.1, ...
                 'Comp_i' , [0, 0, 1], ...
                 'cp'     , 0, ...
                 'name'   , 'PROD', ...
                 'Sign'   , -1);


%% Setup the Schedule
%
% We simulate the formation of a polymer plug.
% Three periods:
% 1) water only         (100 days)
% 2) water + polymer    (50 days)
% 3) water only         (150 days)
%
% We set up a short schedule so that the computations do not take to much
% time in this example.
%

control(1).W = W;
[W([W.sign] > 0).cp] = 2*kilogram/meter^3;
control(2).W = W;

polyinj_stop_time  = 150*day;
end_time           = 300*day;

dt = 10*day;

val1 = linspace(0, polyinj_stop_time, round(polyinj_stop_time/dt));
val2 = linspace(polyinj_stop_time, end_time, round((end_time -polyinj_stop_time)/dt));
step.val     = [diff(val1'); ...
                diff(val2')];
step.control = [2*ones(numel(val1)-1, 1); ... 
                1*ones(numel(val2)-1, 1)];
schedule.step    = step;
schedule.control = control;

% Refine schedule at start.
schedule = refineSchedule(0, day*ones(10, 1), schedule);

%% Setup the initial state
%
nc = G.cells.num;
state0.pressure = ones(nc, 1)*bhpProd;
state0.s        = ones(nc, 1)*[0, 0, 1];
state0.rs       = 0.5*fluid.rsSat(state0.pressure);
state0.rv       = zeros(nc, 1);
state0.cp       = zeros(G.cells.num, 1);
state0.cpmax    = state0.cp;

%% visualize the model properties
%
example_name = 'spe10';
vizPolymerModel();

%% Run the simulation
%
% Generates function handle that set up dynamic plotting for the simulation
fn = getPlotAfterStep(state0, model, schedule, 'plotWell', true, ...
                      'plotReservoir', false);
[wellSols, states] = simulateScheduleAD(state0, model, schedule, 'afterStepFn', ...
                                        fn);


%% Inspect the results
%
figure, plotToolbar(G, states,'field','s:1'); plotWell(G,W,'height',.5);
view(-10,40); axis tight

%% Copyright notice

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
