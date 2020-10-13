mrstModule add dg ad-props ad-core ad-blackoil blackoil-sequential ...
    mrst-gui matlab_bgl vemmech vista
mrstVerbose on

%%

setup = getDGTestCase('simple1d', 'n', 20);

%%

[wsFV, stFV, repFV] = simulateScheduleAD(setup.state0, setup.modelFV, setup.schedule);

%%

[wsDG, stDG, repDG] = deal(cell(numel(setup.modelDG),1));
for dNo = 1:numel(setup.modelDG)
    [wsDG{dNo}, stDG{dNo}, repDG{dNo}] = simulateScheduleAD(setup.state0, setup.modelDG{dNo}, setup.schedule);
end

%%

close all

nc = setup.modelFV.G.cells.num;
x = linspace(1/nc, 1 - 1/nc, nc);

ix = 50;
hold on
plot(x, stFV{ix}.s(:,1));
marker = {'o', 'sq', '^'};
for dNo = 1:numel(setup.modelDG)
    plot(x, stDG{dNo}{ix}.s(:,1), marker{dNo});
end
hold off

%%

close all
dNo = 3;
plotToolbar(setup.modelFV.G, stDG{dNo}, 'plot1d', true);

%% Copyright Notice
%
% <html>
% <p><font size="-1">
% Copyright 2009-2020 SINTEF Digital, Mathematics & Cybernetics.
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
