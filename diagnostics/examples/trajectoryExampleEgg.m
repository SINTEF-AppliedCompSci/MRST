%% trajectoryExampleEgg
% This example launches the (limited robustness) TrajectoryGUI for the Egg 
% model. The GUI enables drawing and draging around well-trajectories which 
% can be saves as an MRST well structues. Diagnostics for the newly created 
% wells can be investigated by launcing the DiagnosticsViewer. 
% 
% INSTRUCTIONS
%
% To create a new well:
% * press "New"
% * right-click inside left-side axis for drawing options ("along line" is recommend). 
% * This enables you to draw the horizontal path of the well. Use left-click on the points 
%   to move them, add additional points by left-clicking where you want new points, and 
%   right-click on the points to delete.
% * On the right-hand axis you now see a vertical cross section along the horizontal path 
%   from the left-hand axis. To draw the path of the well onto this cross section , 
%   right-click in the right-hand axis and choose "Draw points -> along line".
%
% To view the well(s):
% * When you have drawn both a horizontal and cross-sectional well path, you can press the 
%   "View" button to get a 3D view of the well(s) with the traversered grid cells.
%
% Configure and save well:
% * Specify our well parameters by selecting Type, Control and Value
% * Save your new well by pressing the "Save" button. It now appears in the topmost dropdown list
% 
% Launch diagnostics:
% * When you are happy with your well(s) and have saved them, launch diagnostics by pressing the 
%   "Launch diagnostics" button. 

mrstModule add diagnostics mrst-gui ad-core ad-blackoil wellpaths example-suite

gravity off
tmp = setupEggFn(1);
[W, model, state0] = deal(tmp.W, tmp.model, tmp.state0);
model.G = computeGeometry(model.G);

d = TrajectoryGUI(model, W, 'state0', state0);

%%
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
