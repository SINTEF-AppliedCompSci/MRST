mrstModule add dg ad-props ad-core ad-blackoil blackoil-sequential ...
    mrst-gui matlab_bgl vem vemmech vista upr
mrstVerbose on

%%

setup111 = getDGTestCase('qfs_wog_3d', 'useUnstructCubature', true, 'degree', [1,1,1], 'n', 10);
setup110 = getDGTestCase('qfs_wog_3d', 'useUnstructCubature', true, 'degree', [1,1,0], 'n', 10);

%%

[wsFV, stFV, repFV] = simulateScheduleAD(setup111.state0, setup111.modelFV{1}, setup111.schedule);

%%

[wsDG111, stDG111, rep111] = simulateScheduleAD(setup111.state0, setup111.modelDG{1}, setup111.schedule);

%%

[wsDG110, stDG110, rep110] = simulateScheduleAD(setup110.state0, setup110.modelDG{1}, setup110.schedule);

%%

plotWellSols({wsDG111, wsDG111})

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
