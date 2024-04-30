% Example same as blackoilTutorialGravSeg in MRST simulating gravity
% segregation.

clear all
close all

clearSimulation = false;
mrstVerbose off;

mrstModule add ad-core ad-blackoil ad-props mrst-gui mpfa nfvm

meshdiscfactor = 2
timediscfactor = 2

input.twist  = true;
input.useAMG = true;

% Grid
input.dims     = [4, 4, 10] * meshdiscfactor;
input.physDims = [2, 2, 5] * meter;
G = cartGrid(input.dims, input.physDims);

if input.twist
    G = twister(G);
end
G = computeGeometry(G);

% Rock
rock = makeRock(G, 100*milli*darcy, 0.5);
pv = sum(poreVolume(G, rock));

% Fluid
fluid = initSimpleADIFluid('phases', 'WO', ...
                           'mu', [1, 10]*centi*poise, ...
                           'n',  [1, 1], ...
                           'rho', [1000, 700]*kilogram/meter^3);

% Initial state
lower = G.cells.centroids(:, 3) > 0.5*input.physDims(3);
sW = ones(G.cells.num, 1);
sW(lower) = 0;
s = [sW, 1 - sW];
state0 = initResSol(G, 100*barsa, s);

% Plot
figs = 1:5;
plotter(G, state0, 'Initial', figs(1), 0.2);

% BCs
bc = [];
bc = pside(bc, G, 'ZMin', 100*barsa, 'sat', [0 1]);

% Schedule
input.endtime = 750*day;
input.timesteps = 20*timediscfactor;
dt = repmat(input.endtime/input.timesteps, input.timesteps, 1);
schedule = simpleSchedule(dt, 'bc', bc);

% Model
gravity reset on;
model = TwoPhaseOilWaterModel(G, rock, fluid);

% Run
opts = {'clearSimulation', clearSimulation};

input.discmethod = 'TPFA';
statesTPFA = solve(state0, model, schedule, input, opts{:});
plotter(G, statesTPFA{end}, sprintf('%s', input.discmethod), figs(2));

input.discmethod = 'AvgMPFA';
model_avgmpfa = setAvgMPFADiscretization(model, 'myRatio', []);
statesAvgMPFA = solve(state0, model_avgmpfa, schedule, input, opts{:});
plotter(G, statesAvgMPFA{end}, sprintf('%s', input.discmethod), figs(3))

input.discmethod = 'MPFA';
model_mpfa = setMPFADiscretization(model);
statesMPFA = solve(state0, model_mpfa, schedule, input, opts{:});
plotter(G, statesMPFA{end}, sprintf('%s', input.discmethod), figs(4))

input.discmethod = 'NTPFA';
model_ntpfa = setNTPFADiscretization(model, 'myRatio', []);
statesNTPFA = solve(state0, model_ntpfa, schedule, input, opts{:});
plotter(G, statesNTPFA{end}, sprintf('%s', input.discmethod), figs(5))

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

%% Helper Functions

function plotter(G, state, str, figno, edgealpha)

    if nargin == 4
        edgealpha = 0.0;
    end

    subplot(1, 5, figno)
    plotCellData(G, state.s(:,1), 'edgealpha', edgealpha);
    axis equal tight off
    view(50, 20)
    clim([0, 1])
    title(str)
    drawnow

    if figno == 5
        x = get(subplot(1, 5, figno), 'Position');
        colorbar('Position', [x(1)+0.15 x(2)+0.24 0.025 x(4)-0.45])
    end

end


function states = solve(initstate, model, schedule, input, varargin)

    opt = struct('clearSimulation', false, ...
                 'basename', 'testgravity2');
    opt = merge_options(opt, varargin{:});

    if input.useAMG
        threshold = 1;
    else
        threshold = inf;
    end

    ls  = selectLinearSolverAD(model, 'BackslashThreshold', threshold);
    nls = NonLinearSolver('LinearSolver', ls);

    name      = md5sum(input);
    directory = fullfile(mrstOutputDirectory(), opt.basename);

    problem = packSimulationProblem(initstate, model, schedule, opt.basename, ...
                                    'Directory', directory, ...
                                    'Name', name, ...
                                    'NonLinearSolver', nls);

    filename = fullfile(directory, problem.Name, 'input.mat');
    save(filename, 'input');

    if opt.clearSimulation
        clearPackedSimulatorOutput(problem, 'Prompt', false);
    end

    simulatePackedProblem(problem);

    if nargout > 0
        [~, states] = getPackedSimulatorOutput(problem);
    end

end
