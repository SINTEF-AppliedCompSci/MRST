%% Workflow example for MRST-AD
% This example aims to show the complete workflow for creating, running and
% analyzing a simulation model. Unlike the other examples, we will create
% all features of the model manually to get a self-contained script without
% any input files required.
%
% The model we setup is a slightly compressible two-phase oil/water model
% with multiple wells. The reservoir has a layered strategraphy and
% contains four intersecting faults.
%
% Note that this example features a simple conceptual model designed to
% show the workflow rather than a problem representing a realistic scenario
% in terms of well locations and fluid physics.
mrstModule add ad-core ad-blackoil ad-props mrst-gui
close all;

%% Reservoir geometry and petrophysical properties
% We begin by setting up the grid and rock structures. The grid is created
% by "makeModel3", which creates a structured model with intersecting
% faults. We assume a layered permeability structure with 300, 100, and 500
% md in the lower, middle, and top layers. respectively.

% Define grid
% grdecl = makeModel3([50, 50, 5], [1000, 1000, 5]*meter);
% G = processGRDECL(grdecl);
if 0
    G = cartGrid([100, 100, 1], [1000, 1000, 30]*meter);
    G = computeGeometry(G);

    % Set up permeability based on K-indices
    [I, J, K] = gridLogicalIndices(G);

    top = K < G.cartDims(3)/3;
    lower = K > 2*G.cartDims(3)/3;
    middle = ~(lower | top);

    px = ones(G.cells.num, 1);
    px(lower) = 300*milli*darcy;
    px(middle) = 100*milli*darcy;
    px(top) = 500*milli*darcy;

    % Introduce anisotropy by setting K_x = 10*K_z.
    perm = [px, px, 0.1*px];
    rock = makeRock(G, perm, 0.3);
else
    mrstModule add spe10
    [G, ~, rock] = getSPE10setup(1);
    rock.poro(rock.poro < 0.01) = 0.01;
end
figure(1); clf
plotCellData(G, log10(rock.perm(:, 1)), 'edgecolor', 'none')
axis equal tight
%% Define wells and simulation schedule
simTime = 10*year;
nstep   = 50;
refine  = 5;

% Producers
pv      = poreVolume(G, rock);
injRate = 1*sum(pv)/simTime;

offset  = 0;

W = verticalWell([], G, rock, 1, 1, [],...
                'Name', 'P1', 'comp_i', [1 0], 'Val', 250*barsa, 'Type', 'bhp');            
W = verticalWell(W, G, rock,  1,  G.cartDims(2), [],...
                'Name', 'P2', 'comp_i', [1 0], 'Val', 250*barsa, 'Type', 'bhp');
W = verticalWell(W, G, rock, G.cartDims(1), 1, [], ...
                'Name', 'P3', 'comp_i', [1 0], 'Val', 250*barsa, 'Type', 'bhp');
W = verticalWell(W, G, rock, G.cartDims(1), G.cartDims(2), [], ...
                'Name', 'P4', 'comp_i', [1 0], 'Val', 250*barsa, 'Type', 'bhp');

% Injectors
W = verticalWell(W, G, rock, ceil(G.cartDims(1)/2), ceil(G.cartDims(2)/2), 1,...
                'Name', 'I1', 'comp_i', [1 0], 'Val', injRate, 'Type', 'rate');

% W = verticalWell([], G, rock, offset, offset, [],...
%                 'Name', 'P1', 'comp_i', [1 0], 'Val', 250*barsa, 'Type', 'bhp');            
% W = verticalWell(W, G, rock,  offset, floor(G.cartDims(1)/2)+3, [],...
%                 'Name', 'P2', 'comp_i', [1 0], 'Val', 250*barsa, 'Type', 'bhp');
% W = verticalWell(W, G, rock, offset, G.cartDims(2) - offset/2, [], ...
%                 'Name', 'P3', 'comp_i', [1 0], 'Val', 250*barsa, 'Type', 'bhp');
% % Injectors
% W = verticalWell(W, G, rock, G.cartDims(1)-5, offset, 1,...
%                 'Name', 'I1', 'comp_i', [1 0], 'Val', injRate, 'Type', 'rate');

% Compute the timesteps
% startSteps = repmat((simTime/(nstep + 1))/refine, refine, 1);
% restSteps =  repmat(simTime/(nstep + 1), nstep, 1);
% timesteps = [startSteps; restSteps];
timesteps = rampupTimesteps(simTime, simTime/nstep);
% Set up the schedule containing both the wells and the timesteps
schedule = simpleSchedule(timesteps, 'W', W);


%% Set up simulation model

% Three-phase template model
fluid = initSimpleADIFluid('mu',    [1, 1, 0]*centi*poise, ...
                           'rho',   [1000, 700, 0]*kilogram/meter^3, ...
                           'n',     [2, 3, 0]);
if 0
    % Constant oil compressibility
    c        = 0.001/barsa;
    p_ref    = 300*barsa;
    fluid.bO = @(p) exp((p - p_ref)*c);
    clf
    p0 = (100:10:500)*barsa;
    plot(p0/barsa, fluid.bO(p0))
    xlabel('Pressure (bar)')
    ylabel('Ratio')
    title('Reciprocal formation volume factor for oil (bO)')
end

% Construct reservoir model
gravity reset on
model = TwoPhaseOilWaterModel(G, rock, fluid);

%% Define initial state
sW = zeros(G.cells.num, 1);
sat = [sW, 1 - sW];

g = model.gravity(3);
p_res = 300*barsa;
state0 = initResSol(G, p_res, sat);

%% Simulate base case
[wellSols, states, report] = ...
   simulateScheduleAD(state0, model, schedule);



%% Create an upscaled, coarser model
mrstModule add coarsegrid
% cdims = [10, 10, 1];
cdims = ceil(G.cartDims./[10, 10, 1]);
p = partitionUI(G, cdims);
CG = generateCoarseGrid(G, p);

%% Upscale the model and run the coarser problem
mrstModule add incomp agglom upscaling

%%
model_c0 = upscaleModelTPFA(model, p);
state0_c = upscaleState(model_c0, model, state0);
schedule_c0 = upscaleSchedule(model_c0, schedule);
[wellSols_c0, states_c0, rep_c0] = simulateScheduleAD(state0_c, model_c0, schedule_c0);

%%
mrstModule add upscaling relperm-upscale


if 0
    mrstModule add sequential
    f2 = initSimpleADIFluid('mu', [1, 1, 1]*centi*poise);
    pmodel = PressureOilWaterModel(G, rock, f2);
    [hf, T_c, W_c, rep2] = upscaleTrans(CG, pmodel, ...
        'Wells', W, 'bc_method', 'wells', 'fix_trans', true, 'state', state0, 'dt', 1*day, 'match_method', 'max_flux');
    model_c = upscaleModelTPFA(model, p, 'transCoarse', T_c);


    schedule_c = schedule;
    for i = 1:numel(schedule_c.control)
        schedule_c.control(i).W = W_c;
    end
else
    ts = computeTransmissibilityFromStates(p, states, model, schedule);
    [model_c, schedule_c] = assignCoarseTransmissibility(p, model, schedule, ts);
end

[wellSols_c, states_c, rep_c] = simulateScheduleAD(state0_c, model_c, schedule_c);
%%
kr0 = computeRelpermFromStates(states, model_c, model, schedule_c, schedule);

kr = kr0;
[model_kr, kr] = assignUpscaledRelperm(model_c, kr, 'extrapolateend', true, 'monotoneMethod', 'remove', 'setWells', false);

[wellSols_kr, states_kr, rep_kr] = simulateScheduleAD(state0_c, model_kr, schedule_c);
%%
kr = kr0;
% kr = mergeHalfFaceRelperm(model_c, kr0, 'type', 'face');
[model_krf, krf] = assignUpscaledRelperm(model_c, kr, 'extrapolateend', true, 'monotoneMethod', 'slope', 'setWells', true);

[wellSols_krf, states_krf, rep_krf] = simulateScheduleAD(state0_c, model_krf, schedule_c);

%%
plotWellSols({wellSols, wellSols_c0, wellSols_c, wellSols_kr, wellSols_krf}, report.ReservoirTime, 'datasetnames', {'Fine', 'Harmonic', 'T upscale', 'T + kr upscale', 'T + kr upscale (with wells)'})
%%
%%
kr = kr0
close all
for i = 1:2
    figure;
    hold on
    s = 0:0.01:1;
    if i == 1
        k = model.fluid.krW(s);
    else
        k = model.fluid.krO(s);
    end
    
    
    for j = 1:numel(kr{i}.reservoir)
        d = kr{i}.reservoir{j};
        
        xx = d.S;
        ff = d.kr;
        if 1
            [xx, ff] = regularizeSaturationFunction(xx, ff, 'monotone', true, 'extrapolateEnd', true, 'monotoneMethod', 'slope');
        end
        if any(d.S)
            plot(xx, ff, '-o')
%             j
%             pause()
        end
    end
    plot(s, k, 'b', 'linewidth', 2);
end

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
