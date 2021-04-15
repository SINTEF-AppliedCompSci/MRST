mrstModule add dg vem vemmech ad-props ad-core ad-blackoil ...
    blackoil-sequential mrst-gui reorder matlab_bgl ...
    ad-eor trust-region weno deckformat
mrstVerbose on
gravity reset off

%% Set up fluid

if 1
    fluid = initSimpleADIFluid('phases', 'WO'                      , ...
                               'rho'   , [1000, 1]*kilogram/meter^3, ...
                               'mu'    , [0.5 , 1]*centi*poise     , ...
                               'n'     , [1   , 1]                 );
    fluid = addSimplePolymerProperties(fluid, 'cmax', 5);
else
    fn   = fullfile(mrstPath('dg'), 'examples', 'siamgs-19', 'data', 'POLYMERDATA.DATA');
    deck = readEclipseDeck(fn);
    %The deck variables are converted to SI units.
    deck = convertDeckUnits(deck);
    % Setup the fluid structure
    fluid = initDeckADIFluid(deck);
    fluid.effads = @(c, cmax) effads(c, cmax, fluid);
end


%% Set up schedule

% Time and rate
l     = 1000*meter;
w     = 0.01*l;
poro  = 0.4;
time  = 0.5*year;
wTime1 = 0.5*year;
wTime2 = 5*year;
pTime = 0.3*year;
rate  = 0.2*l*w*poro/year;
dt    = 10*day;


%%

dataDir  = fullfile(mrstPath('dg'), 'examples', 'siamgs-19', 'data', 'polymer1d');
pack = @(name, n, state0, model, schedule) ...
            packSimulationProblem(state0, model, schedule, num2str(n), ...
                           'Directory', fullfile(dataDir, num2str(n)), ...
                           'Name'     , name                         );

%%

n      = [5, 10, 50, 100];
nG     = numel(n);
degree = 0:4;
[problemFV, problemWENO] = deal(cell(nG,1));
problemDG = cell(nG,numel(degree));
for gNo = 1:nG
    G     = computeGeometry(cartGrid([n(gNo),1], [1, 0.01]*l));
    G     = computeCellDimensions2(G);
    rock  = makeRock(G, 100*milli*darcy, 0.4);
    % Wells
    W = [];
    W = addWell(W, G, rock, 1   , ...
                'type'  , 'rate', ...
                'val'   , rate  , ...
                'comp_i', [1,0] );
    W = addWell(W, G, rock, G.cells.num, ...
                'type'  , 'bhp'   , ...
                'val'   , 50*barsa, ...
                'comp_i', [1,0]   );
    [W.c]    = deal(0);
    schedule = makePolymerSlugSchedule(W, fluid, 'dt', 10*day, 'pTime', pTime, 'wTime', [wTime1, wTime2]);
    % Initial state
    sW          = 0.0;
    state0      = initResSol(G, 100*barsa, [sW,1-sW]);
    [state0.c, state0.cmax] = deal(zeros([G.cells.num, 1]));
    % FI model
    model = OilWaterPolymerModel(G, rock, fluid);
    % Sequential problem
    modelFV = getSequentialModelFromFI(model);
    problemFV{gNo} = pack('fv', n(gNo), state0, modelFV, schedule);
    % WENO problem
    if 0
        fluid.effads  = @(c, cmax) effads(c, cmax, fluid);
        fluid.krOW    = fluid.krO;
        fluid.pvMultR = @(p) 1;    
        modelWENO = PolymerHighOrderModel(G, rock, fluid, []);
        modelWENO.FacilityModel %= PolymerHighOrderModelWell(G, rock, fluid, []);
    else
        modelWENO = HigherOrderOilWaterPolymerModel(G, rock, fluid, []);
        weno      = WENODiscretization(model, G.griddim, ...
                      'includeBoundary'     , true, ...
                      'interpolateReference', true);
        modelWENO.FluxDiscretization.saturationDiscretization = weno;
        modelWENO.FluxDiscretization.relPermDiscretization = weno;
        modelWENO.FluxDiscretization.discritizeRelPerm = false;
        modelWENO.FluxDiscretization.discritizeViscosity = false;
        problemWENO{gNo} = pack('weno', n(gNo), state0, modelWENO, schedule);
    end

    % DG problems
    for k = degree
        modelDG = modelFV;
        modelDG.transportModel ...
            = TransportOilWaterPolymerModelDG(G, rock, fluid, ...
                                       'dsMaxAbs'     , 0.2 , ...
                                       'degree'       , k   , ...
                                       'meanTolerance', 1e-3, ...
                                       'outTolerance' , 1e-3);
        state0 = assignDofFromState(modelDG.transportModel.disc, state0);
        problemDG{gNo, k+1} = pack(['dg-', num2str(k)], n(gNo), state0, modelDG, schedule);
    end
    
    
end

%%

nRef  = 10000;
G     = computeGeometry(cartGrid([nRef,1], [1, 0.01]*l));
rock  = makeRock(G, 100*milli*darcy, 0.4);
% Wells
W = [];
W = addWell(W, G, rock, 1   , ...
            'type'  , 'rate', ...
            'val'   , rate  , ...
            'comp_i', [1,0] );
W = addWell(W, G, rock, G.cells.num, ...
            'type'  , 'bhp'   , ...
            'val'   , 50*barsa, ...
            'comp_i', [1,0]   );
[W.c]    = deal(0);
schedule = makePolymerSlugSchedule(W, fluid, 'dt', 10*day, 'pTime', pTime, 'wTime', [wTime1,wTime2]);
% Initial state
sW          = 0.0;
state0      = initResSol(G, 100*barsa, [sW,1-sW]);
[state0.c, state0.cmax] = deal(zeros([G.cells.num, 1]));
% FI model
model = OilWaterPolymerModel(G, rock, fluid);
% Sequential problem
problemRef = pack('fi', nRef, state0, model, schedule);

%% Simulate reference problem

simulatePackedProblem(problemRef);

%% Simulate FV problems

simulatePackedProblem(problemFV);

%% Simulate WENO problems

simulatePackedProblem(problemWENO);

%% Simulate DG problems

simulatePackedProblem(problemDG);

%% Debug

% setup = problemDG{1,2}.SimulatorSetup;
setup = problemWENO{1}.SimulatorSetup;
% setup = problemFV{1}.SimulatorSetup;
% setup = struct('model', model, 'state0', state0, 'schedule', schedule);
[ws, st, rep] = simulateScheduleAD(setup.state0, setup.model, setup.schedule);

%% Load reference results

[wsRef, stRef, repRef] = getPackedSimulatorOutput(problemRef);

%% Load FV results

[wsFV, stFV, repFV] = getMultiplePackedSimulatorOutputs(problemFV);

%%

[wsWENO, stWENO, repWENO] = getMultiplePackedSimulatorOutputs(problemWENO);

%% Load dG results

[wsDG, stDG, repDG] = deal(cell(nG,numel(degree)));
for k = degree
    [wsDG(:,k+1), stDG(:,k+1), repDG(:,k+1)] = getMultiplePackedSimulatorOutputs(problemDG(:,k+1));
end

%%

get_wcut = @(ws) cellfun(@(w) w(2).wcut, ws);
wcutRef = get_wcut(wsRef);
clr    = lines(5);
refClr = clr(2,:);
clr    = clr([1, 3:end],:);

close all
tvec = cumsum(schedule.step.val)/day;
for gNo = 1:nG
    figure
    hold on
    wcutDG = zeros(numel(tvec), numel(degree));
    for dNo = 1:2
        wcutDG(:, dNo) = get_wcut(wsDG{gNo, dNo});
        plot(tvec, wcutDG(:, dNo), 'color', clr(dNo,:), 'linew', 1)
    end
    wcutWENO = get_wcut(wsWENO{gNo});
    plot(tvec(1:numel(wcutWENO)), wcutWENO, '.', 'color', clr(:,1))
    plot(tvec, wcutRef, '--', 'color', refClr, 'linew', 4)
    d = 0.01;
    axis([0, tvec(end), -d, 1+d])
    box on
    hold off
end

%%

close all

gNo = 4;
k   = 1;

% setup = problemDG{gNo,k+1}.SimulatorSetup;
% st    = stDG{gNo, k+1};

setup = problemFV{gNo}.SimulatorSetup;
st    = stFV{gNo};

% setup = problemRef.SimulatorSetup;
% st    = stRef;

plotToolbar(setup.model.G, st, 'plot1d', true);

%%

close all
plotWellSols([wsFV, wsDG(:,2)], schedule.step.val)

%%

if 1
    clearPackedSimulatorOutput(problemRef, 'prompt', false);
end

%%

if 1
    clearPackedSimulatorOutput(problemFV, 'prompt', false);
end

%%

if 1
    clearPackedSimulatorOutput(problemDG, 'prompt', false);
end

%% Copyright Notice
%
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
 
