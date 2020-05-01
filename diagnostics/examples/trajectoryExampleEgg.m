%% trajectoryExampleEgg
% This example launches the (limited robustness) TrajectoryGUI for the Egg 
% model. The GUI enables drawing and draging around well-trajectories which 
% can be saves as an MRST well structues. Diagnostics for the newly created 
% wells can be investigated by launcing the DiagnosticsViewer. 
% 
% INSTRUCTIONS:
% *  press new to create new well
% *  right-click inside axis for drawing options ("along line" is recommend)
% *  left-click edge of points to drag
% *  right-click (edge of) points for options
% *  click view to plot resulting trajectory in 3D and traversed grid cells
% *  click save to save well(s)
% *  click launch diagnsotics to view diagnostics for saved well

mrstModule add diagnostics mrst-gui ad-core ad-blackoil wellpaths

gravity off
tmp = setupEggFn(1);
[W, model] = deal(tmp.W, tmp.model);
model.G = computeGeometry(model.G);

d = TrajectoryGUI(model, W);

%%
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
