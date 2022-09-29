%% trajectoryExampleEgg
% This example launches the TrajectoryGUI for the Egg model. The GUI enables 
% editing, drawing and draging around well-trajectories which can be saved as 
% MRST well structues. Diagnostics for the newly edited/created wells can 
% be investigated by launcing the DiagnosticsViewer. 
% 
% INSTRUCTIONS for editing a well-trajectory and comparing diagnostics with
% the base-case.
% 1) Make a new well-case by clicking "New" in the "Well cases"-box. We now 
%    have two identical well-cases ("base" and "case 1").
% 2) In the "Edit/Add well"-box, select a well-name in the drop-down menu.
% 3) Edit the well-trajectory by 
%    a) Drag the heel/toe green endpoints along the 2D-slice in the left plot. 
%       Change the slice by adjusting the slice-angle in the "Trajectory slice
%       options"-box.
%    b) To make curved trajectories, add points along the trajectory:
%       Right-click inside axis of the left plot -> Draw points -> Along line. 
%       Left-click to place new points. Right-click point to delete. To switch 
%       off point-drawing: Right-click inside axis of the left plot -> Draw 
%       points -> off.
%    c) To move the entire trajectory in the horizontal plane, drag pink 
%       square in the right plot.
%    d) To undo all changes, press "Undo".
% 4) Save the new well-trajectory in "case 1" by pressing "Save".
% 5) Optionally modify and save other wells in "case 1".
% 6) Optionally make additional cases by repeating steps 1)-5)
% 7) Select all cases listed in the "Well cases"-box, and press 
%    "Launch diagnostics" to compare e.g., the effect an updated trajectory 
%    has on sweep- and drainage-regions
%
% Other options in GUI, includes drawing new wells, and changing well control 
% types and values (the latter may produce strange results unless values are
% selected with care). 

mrstModule add diagnostics mrst-gui ad-core ad-blackoil wellpaths test-suite

gravity off
tmp = setupEggFn(1);
[W, model, state0] = deal(tmp.W, tmp.model, tmp.state0);
model.G = computeGeometry(model.G);

d = TrajectoryGUI(model, W, 'state0', state0);

%%
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
