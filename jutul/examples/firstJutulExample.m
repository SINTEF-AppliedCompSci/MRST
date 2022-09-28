%% First introductory example to Jutul as a MRST accelerator
% JutulDarcy is a reservoir simulator written in Julia by SINTEF Digital
% based on the Jutul solver framework, developed for high-performance
% testing of new simulator concepts. We can thus set up a case in MRST and
% run it in Jutul for increased computational performance. Jutul supports
% black-oil, immiscible, and compositional models with multisegment wells.
%
% To install Jutul, you must first install the latest version of Julia [1].
% Once installed, run Julia and add the JutulDarcy package. If you are
% interested in using Julia for other things, we recommend adding it to an
% environment[2]. Otherwise, you can add it to the default environment by
% running the following command in the Julia prompt:
%
% using Pkg; Pkg.add("JutulDarcy")
%
% Once downloaded, you are ready to run this example. Note that the example
% pauses once the simulation is ready to be run in the Julia terminal.
%
% For more details on Jutul, see the JutulDarcy repository on GitHub [3]

% [1] https://julialang.org/downloads/
% [2] https://pkgdocs.julialang.org/v1/environments/
% [3] https://github.com/sintefmath/JutulDarcy.jl

mrstModule add ad-core ad-blackoil spe10 deckformat ad-props test-suite compositional jutul
if ~exist('name', 'var')
    name = 'qfs_wo';
end
%%
switch name
    case 'qfs_wo'
        setup = qfs_wo();
    case 'spe10_layer'
        setup = spe10_wo('layers', 1);
    case 'fractures_compositional'
        setup = fractures_compositional();
        % Jutul has built-in separator support. We set up a single
        % separator with conditions that match the default in the other
        % simulator.
        setup.model.AutoDiffBackend = AutoDiffBackend();
        s  = EOSSeparator('pressure', 101325.0, 'T', 288.15);
        sg = SeparatorGroup(s);
        sg.mode = 'moles';
        setup.model.FacilityModel.SeparatorGroup = sg;
    otherwise
        % Assume some test-suite case
        setup = eval(name);
end
[state0, model, schedule] = deal(setup.state0, setup.model, setup.schedule);

%% Write case to disk and output Jutul command
% The routine will write the MRST setup to disk in a *.mat file and then
% output a command you can paste into an existing Julia session that
% already has JutulDarcy preloaded. Once the resulting simulation is
% finished, you can hit the return button to continue execution in MATLAB.
% Notice that the first time you run the simulation, it will take a long
% time since Julia has to compile necessary code. Once compiled, however,
% the simulation is fast.
pth = writeJutulInput(state0, model, schedule, name);
disp('Pausing - run the command in Julia and hit any key to continue')
pause()

%% Once simulated, read back as MRST format
[ws, states] = readJutulOutput(pth);

%% Simulate MRST for comparison purposes
nls = getNonLinearSolver(model);
[ws_m, states_m] = simulateScheduleAD(state0, model, schedule, 'NonLinearSolver', nls);

%% Compare the results
% The two simulators use somewhat different well models: MRST is set up
% with a standard Peacemann well model while Jutul uses multisegment wells
% by default based on tables generated from MRST.  This can give certain
% differences for high flow rates, e.g., in the near-well region. This is
% particularly evident during the initial transitional phase for the SPE10
% example.
mrstModule add mrst-gui
G = model.G;
figure;
plotToolbar(G, states_m, 'field', 'pressure')
title('MRST')
figure;
plotToolbar(G, states, 'field', 'pressure')
title('Jutul')
figure;
plotToolbar(G, applyFunction(@(x, y) compareStructs(x, y), states, states_m), 'field','pressure')
title('Difference')

%% Plot and compare wells
% Results will differ because the two simulators use different well models
plotWellSols({ws, ws_m}, cumsum(schedule.step.val), 'datasetnames', {'Jutul', 'MRST'})

%% Copyright Notice
%
% <html>
% <p><font size="-1">
% Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.
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
