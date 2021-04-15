function W = addSleipnerWellsTrajectory(G,rock,rate)
% Function to create injection well for Sleipner 2019 benchmark model
%
% SYNOPSIS:
%   W = addSleipnerWellsTrajectory(G,rock,rate);
%
%
% PARAMETERS:
%  G        - Finescale Grid of Multilayered Sleipner Benchmark model without
%              caprock cells. Output from getMultiLayerSleipnerGrid.
%
%  rock     - rock structure corresponding to G.
%
%  rate     - Injection rate.  
%
% RETURNS:
%  feederCells - Well structure for well 15_9_A16 from Sleipner 2019 
%                   benchmark model
%
%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}
%%
require wellpaths
% Prior to accurate intersecting algoritms, a more cost effective rough search 
% is performed based on grid face bounding boxes. To avoid overhead, add these 
% fields to the grid prior to calling 'computeTraversedCells'
G = addBoundingBoxFields(G);


%% Define well trajectory and compute intersection with grid


datafolder = fullfile(getDatasetPath('sleipner2019'),...
    'well_data','data','Wells_released_2011');

% Get well trajectories

wellDataFile = fullfile(datafolder,'Well 159_A16','15_9_A16_wellposition');
wdata = dlmread(wellDataFile,'',24,0);
pnts = wdata(:,[1 2 7]);
pnts(:,3) = pnts(:,3).*-1;

perfDepth = [1010.5,1013.3];
traj = interp1(pnts(:,3),pnts,(0:.01:1).*(perfDepth(2)-perfDepth(1)) + perfDepth(1));


% Compute intersections between grid and trajectory. Returned struct
% contains cell indiices of traversed cells, vector representing part of 
% segment (exit point minus entry point) and a weight which which is typically 
% one except in cases where a segment is shared between multiple cells   
T = computeTraversedCells(G, traj);  %#ok

%% Add wells from trajectory
% Wells can be added directly from trajectory similarly to the addWell-function. 
% The well indices are computed by the formula
%    WI = sqrt( sum_{k=x,y,z} (Lk*WIk/Dk)^2  ), 
% where [Dx, Dy, Dz] are the cell dimensions, [Lx, Ly, Lz] is the vector of
% the traversing segment (difference between exit and entry-point) and 
% WIx, WIy and WIz are the well-indices for fully penetrating wells in each 
% coordinate direction

W = [];
W = addWellFromTrajectory(W, G, rock, traj, 'Name', 'w1', 'Type', 'rate','Val',rate,'Name','15/9_A16','comp_i',[0 1]);

