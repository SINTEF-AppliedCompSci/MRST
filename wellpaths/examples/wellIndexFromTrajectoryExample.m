%% Setup a simple layered model
mrstModule add wellpaths

rng(0);
[nx, ny, nz] = deal(50, 50, 15);
G = computeGeometry(processGRDECL(makeModel3([nx ny nz])));
% Prior to accurate intersecting algoritms, a more cost effective rough search 
% is performed based on grid face bounding boxes. To avoid overhead, add these 
% fields to the grid prior to calling 'computeTraversedCells'
G = addBoundingBoxFields(G);

[~, ~, K] = gridLogicalIndices(G);
poroLayers = .2 + .5*rand(nz,1);
poro       = poroLayers(K);
permX      = (100*milli*darcy)*poro.^3./(1-poro.^2);
rock       = makeRock(G, [permX, permX, .1*permX], poro);
figure, plotCellData(G, log(rock.perm(:,1))); view(3)

%% Define well trajectory and compute intersection with grid
pnts = [100 150 -15
        200 250 -1
        400 300  6
        700 600  10];
    
dp     = diff(pnts);
len    = [0; cumsum(sqrt(dot(dp, dp, 2)))];
% interpolate 100 points along splines passing trough pnts
traj   = interp1(len, pnts, (0:.01:1)*len(end), 'spline');

% Compute intersections between grid and trajectory. Returned struct
% contains cell indiices of traversed cells, vector representing part of 
% segment (exit point minus entry point) and a weight which which is typically 
% one except in cases where a segment is shared between multiple cells   
T = computeTraversedCells(G, traj)  %#ok

% plot trajectory, vertical intersection, and cells traversed by trajectory
figure, hold on
col = log(rock.perm(:,1));
poly = computeVerticalGridIntersection(G, traj);
patch('Faces', poly.nodes, 'Vertices', poly.coords3D, ...
      'FaceVertexCData', col(poly.cellIx), ...
      'FaceColor', 'flat', 'EdgeAlpha', .3, 'FaceAlpha', .2)
plot3(traj(:,1), traj(:,2), traj(:,3), '-b','LineWidth', 2);
plotCellData(G, col, T.cell, 'EdgeAlpha', .5, 'FaceAlpha', .5);
axis off, view([-1 2 2]), daspect([1 1 .1]), camproj perspective

%% Add wells from trajectory
% Wells can be added directly from trajectory similarly to the addWell-function. 
% The well indices are computed by the formula
%    WI = sqrt( sum_{k=x,y,z} (Lk*WIk/Dk)^2  ), 
% where [Dx, Dy, Dz] are the cell dimensions, [Lx, Ly, Lz] is the vector of
% the traversing segment (difference between exit and entry-point) and 
% WIx, WIy and WIz are the well-indices for fully penetrating wells in each 
% coordinate direction

W = [];
W = addWellFromTrajectory(W, G, rock, traj, 'Name', 'w1', 'Type', 'bhp');
% Make another version of the well with corrections for exterior faces. The
% correction approximates the ratio r of radial flow from the segment that 
% is blocked by an exterior face in the well cell, and "corrects" the well
% index by a factor 1-r;
W = addWellFromTrajectory(W, G, rock, traj, 'Name', 'w2', 'Type', 'bhp', ...
                         'exteriorFaceCorrection', true);
%plot cells where this correction has an effect:
hasEffect = W(1).WI>W(2).WI;
figure, hold on
plot3(traj(:,1), traj(:,2), traj(:,3), '-b','LineWidth', 2);
plotGrid(G, 'FaceColor', 'none', 'EdgeAlpha', .1)
plotGrid(G, W(1).cells(hasEffect), 'FaceColor', 'r', 'FaceAlpha', .5);
plotGrid(G, W(1).cells(~hasEffect), 'FaceColor', 'g', 'FaceAlpha', .5);
axis off, view([-1 2 2]), daspect([1 1 .1]), camproj perspective

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
