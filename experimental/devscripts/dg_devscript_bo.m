mrstModule add dg ad-props ad-core ad-blackoil blackoil-sequential ...
    mrst-gui matlab_bgl vem vemmech vista
mrstVerbose on

%%

setup = getDGTestCase('spe1');

%%

ix = 1:7;
schdl = setup.schedule;
schdl.step.val = schdl.step.val(ix);
schdl.step.control= schdl.step.control(ix);
sim = @(model, inx) simulateScheduleAD(setup.state0, setup.(model){inx}, schdl);

%%

[wsFV, stFV, repFV] = sim('modelFV', 1);

%%

[wsDG{1}, stDG{1}, repDG{1}] = sim('modelDG', 1);

%%

[wsDG{2}, stDG{2}, repDG{2}] = sim('modelDG', 2);

%%

[wsFI, stFI, repFI] = sim('modelFI', 1);

%%

plotWellSols({wsFI, wsDG{:}})

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
