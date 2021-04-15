mrstModule add ad-core
%% Conceptual example
group = StateFunctionGrouping('mystorage'); % Store in state.mystorage
% Set the numbers
group = group.setStateFunction('A', TutorialNumber('a'));
group = group.setStateFunction('B', TutorialNumber('b'));
group = group.setStateFunction('X', TutorialNumber('x'));
group = group.setStateFunction('Y', TutorialNumber('y'));
% Create products
group = group.setStateFunction('AB', TutorialProduct('A', 'B'));
group = group.setStateFunction('XY', TutorialProduct('X', 'Y'));
% Create the final sum
group = group.setStateFunction('G', TutorialAdd('AB', 'XY'));
disp(group)
%% Plotting
clf;
[h, g] = plotStateFunctionGroupings(group, 'TextArg', {'FontSize', 16, 'horizontalalignment', 'center', 'verticalalignment', 'middle'}, 'includestate', false);
colormap(brighten(lines, 0.5))
printStateFunctionGroupingTikz(g)
%% Show evaluation
clc
state0 = struct('a', 3, 'b', 5, 'x', 7, 'y', 2); % Initial state
state = group.initStateFunctionContainer(state0); % Enable caching
disp(state.mystorage)

disp('First call:')
G = group.get([], state, 'G'); % First call, triggers execution
disp(state.mystorage)
disp('Second call:')
G = group.get([], state, 'G') % Second call, cached
%% Partial caching
clc
state = group.initStateFunctionContainer(state0);
xy = group.get([], state, 'XY'); % First call, triggers execution
G = group.get([], state, 'G'); % Second call, cached
%% Dependency checking
mrstVerbose on; % Print missing dependencies
group.checkDependencies();
%% Add an unmet dependency
group = group.setStateFunction('F', TutorialProduct('X', 'Z'));
ok = group.checkDependencies();
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
