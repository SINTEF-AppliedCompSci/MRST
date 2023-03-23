%% Example
% This script contains an example a well branching. It contains several 
% well-intersections. 

%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2016-2020 Runar Lie Berge. See COPYRIGHT.TXT for details.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}  

%% Set well paths
% voronoi2D store well as arrays in cells. Each row in 
% an array correspond to a corner of the path.

well = {[0.5,0.2; 0.5,0.3;0.47,0.4;0.4,0.5; 0.33,0.6;0.26,0.7], ...
        [0.5,0.3;0.53,0.4;0.58,0.5],...            
        [0.5,0.45;0.5,0.55;0.45,0.65;0.4,0.75;0.38,0.85],...
        [0.5,0.55;0.55,0.65;0.6,0.75;0.62,0.85]};

% We now plot the wells to see how they look like
figure('Position',[480 340 980 420])
subplot(1,2,1)
plotLinePath(well,'color','blue');
axis equal tight
axis ([0,1,0,1])
title('well paths')
drawnow
%% Setting gridding parameters:
% Before we call the gridding fuctions we set some parameters.
gS  = [1/24,1/24]; % The grid size
wGf = 0.25;        % The relative size of well cells compared to gS
fGf = 0.5;         % The relative size of fault cells compared to gS
nRs = 2;           % Number of refinement steps towards the well
mLs = [0.1,0.05]'; % This sets the distance from a well where each 
                   % each refinement step should start.

%% Create Grid
% We can now create the composite Pebi grid:
Gc = compositePebiGrid2D(gS, [1, 1], ...
                       'cellConstraints', well, 'CCFactor',wGf, ...
                       'mlqtMaxLevel', nRs,'mlqtLevelSteps', mLs);

%% Plot composite grid
%And plot it
plotGrid(Gc, 'facecolor','none')
axis equal tight
title('compositePebiGrid2D(...)')
drawnow
% We plot the cells taged as well cells
plotGrid(Gc,Gc.cells.tag,'facecolor','b')

%% Set pebiGrid2D Parameters:
% We now use the other wrapper function to create a PEBI-grid using
% distmesh:
eps = 1/12; % This parameter defines the refinement around the wells. The
            % cell size increase exponentially from the wells with 
            % exp(-(distance from well)/eps);

%% Generate grid
% distmesh will most likely not converge in the maximum number of iterations.
% This is usually not a problem, since the grid most likely is good before
% the convergence requirement is met. 

Gdist = pebiGrid2D(1/24, [1, 1], 'cellConstraints', well, ...
                'CCFactor', 0.5^2, 'CCRefinement', true, ...
                'CCEps',eps);
%% Plot pebiGrid2D
subplot(1,2,2), hold on
plotGrid(Gdist,'facecolor','none')
axis equal tight
title('pebiGrid2D(...)')
% We plot the cells taged as well cells
plotGrid(Gc,Gc.cells.tag,'facecolor','b')
