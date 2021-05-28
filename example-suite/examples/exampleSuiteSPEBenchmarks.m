%% Minimal example showing how to run multiple examples using MRSTExample

%% Add required modules
mrstModule add example-suite
mrstModule add ad-core ad-props ad-blackoil deckformat
mrstModule add mrst-gui
mrstVerbose on

%% Simulate three SPE benchmarks
examples = {'spe1_bo', 'spe3_bo', 'spe9_bo'};
for ex = examples
    example = MRSTExample(ex{1});
    example.plot(); drawnow, pause(1);
    problem = example.getPackedSimulationProblem();
    simulatePackedProblem(problem);
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
