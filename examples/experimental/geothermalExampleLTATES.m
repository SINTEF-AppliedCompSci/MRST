mrstModule add geothermal
mrstModule add upr libgeometry
mrstModule add ad-core ad-props compositional
mrstModule add test-suite spe10
mrstModule add mrst-gui

mrstVerbose on

%%
test = TestCase('ltates_geothermal');

%%
close all
aquifer = test.model.G.cells.layer == 2 | ...
          test.model.G.cells.layer == 3;
test.figure();
plotCellData(test.model.G, log10(test.model.rock.perm),  aquifer, ...
                                                         'edgeAlpha', 0.5);
plotCellData(test.model.G, log10(test.model.rock.perm), ~aquifer, ...
                                    'edgeColor', 'none', 'faceAlpha', 0.2);
test.plotWells();
test.setAxisProperties(gca); camlight; colormap(pink);


%%
tss = StateChangeTimeStepSelector('firstRampupStepRelative', 0.1, ...
                                  'resetOnControlsChanged', true, ...
                                  'maxRelativeAdjustment', 5, ...
                                  'targetProps', {'T'}, ...
                                  'targetChangeAbs', 5*Kelvin);
nls = NonLinearSolver('TimestepSelector', tss);
problem = test.getPackedSimulationProblem(...
            'NonLinearSolver', nls, ...
            'ExtraArguments', {'controlLogicFn', test.options.ctrlFn});
        
%%
simulatePackedProblem(problem, 'restartStep', 1);

%%
[wellSols, states, reports] = getPackedSimulatorOutput(problem);

%%
for i = 1:numel(states)
    states{i}.cells = (test.model.G.cells.layer == 2 | ...
                       test.model.G.cells.layer == 3)*1;
end

%%
test.plot(states);

%%
plotWellSols(wellSols, test.schedule.step.val);

%%
close all
aquifer = test.model.G.cells.layer == 2 | ...
          test.model.G.cells.layer == 3;
test.figure();
plotCellData(test.model.G, states{40}.T,  aquifer & states{40}.T > 315, 'edgeAlpha', 0.5);
plotGrid(test.model.G, aquifer, 'faceColor', 'none', 'edgeAlpha', 0.2);
test.plotWells();
a = 0.9; cmap = summer; cmap = cmap*a + (1-a);
test.setAxisProperties(gca); camlight; colormap(cmap);

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