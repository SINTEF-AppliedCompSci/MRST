%% Example demonstrating how to run black oil cases from .DATA files in Jutul
% Black oil cases can also be run in Jutul, provided that they are
% constructed from a .DATA file.
mrstModule add ad-core ad-blackoil ad-props deckformat jutul
if ~exist('name', 'var')
    name = 'spe1';
end
if ~exist('use_daemon', 'var')
    use_daemon = false;
end
reorder = {};
switch name
    case 'spe1'
        pth = getDatasetPath('spe1');
        deck_path  = fullfile(pth, 'BENCH_SPE1.DATA');
    case 'spe9'
        pth = getDatasetPath('spe9');
        deck_path  = fullfile(pth, 'BENCH_SPE9.DATA');
    case 'egg'
        deck_path = getDeckEGG();
    otherwise
        error('No such case.')
end
[state0, model, schedule, nls] = initEclipseProblemAD(deck_path, 'ReorderStrategy', reorder);
%% Write case to temporary directory and run
% You can set daemon mode to true if set up, see backgroundJutulExample for more details.
[ws, states] = simulateScheduleJutul(state0, model, schedule, 'name', name, 'daemon', use_daemon);
%% Run in MRST
[ws_m, states_m] = simulateScheduleAD(state0, model, schedule, 'NonLinearSolver', nls);
%% Compare wells
plotWellSols({ws, ws_m}, 'datasetnames', {'Jutul', 'MRST'})

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
