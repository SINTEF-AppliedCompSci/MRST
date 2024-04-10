mrstModule add dg ad-props ad-core ad-blackoil blackoil-sequential ...
    mrst-gui matlab_bgl vem vemmech vista
mrstVerbose on

%%

setup    = getDGTestCase('qfs_wo_2d', 'n', 10, 'nkr', 1);
setup.modelDG{2}.transportModel.formulation = 'missingPhase';
sim = @(model, inx) simulateScheduleAD(setup.state0, setup.(model){inx}, setup.schedule);

%%

[wsFV, stFV, repFV] = sim('modelFV', 1);

%%

[wsDG, stDG, repDG] = sim('modelDG', 3);

%%

[wsFI, stFI, repFI] = sim('modelFI', 1);

%%

plotWellSols({wsFI, wsFV, wsDG});

%%

ix = 3;
disc = setup.modelDG{ix}.transportModel.disc;

close all
for t = 1:numel(stDG)
    plotSaturationDG(disc, stDG{t}, 'edgecolor', 'k', 'edgealpha', 0.2);
    view(3);
    zlim([0,1])
    view([100,50]);
    pause(0.1);
end

%%

G = setup.modelFI{1}.G;
plotToolbar(G, stDG)

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
