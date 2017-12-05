%% Example
% This script contains an example of a single well intersected by several
% faults. It uses the two wrapper-functions compositePebiGrid(..) and 
% pebiGrid(..) to create a PEBI-grid conforming to the faults and wells.

%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2016 Runar Lie Berge. See COPYRIGHT.TXT for details.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}  

%% Set well and fault paths
% voronoi2D store well and fault paths as arrays in cells. Each row in 
% an array correspond to a corner of the path.

% We start by creating a well
x     = linspace(0.1, 0.8, 10);
y     = 0.2 + x.^2;
well  = {[x' , y']};

% And four faults crossing the well.
fault   = {[0.3,0.1;0.1,0.3],...        
          [0.3,0.1;0.4,0.7],...
          [0.6,0.2;0.3,0.7],...
          [0.8,0.4;0.5,0.8]};


% We now plot the well and fault to see how it looks like
figure(); hold on
plotLinePath(well,'color','blue');
plotLinePath(fault,'color','red');
axis equal tight
axis equal tight
axis ([0,1,0,1])
title('well & fault paths')
drawnow
axis([0,1,0,1])

%% Setting gridding parameters:
% Before we call the gridding fuctions we set some parameters.
gS  = 1/24; % The grid size
wGf = 0.5;  % The relative size of well cells compared to gS
fGf = 0.5;  % The relative size of fault cells compared to gS
nRs = 1;    % Number of refinement steps towards the well

%% Create Grid
% We can now create the composite Pebi grid:
Gc = compositePebiGrid([gS,gS], [1, 1], ...
                       'wellLines', well, 'wellGridFactor',wGf, ...
                       'faultLines',fault, 'faultGridFactor', fGf,...
                       'mlqtMaxLevel', nRs);

%% Plot composite grid
%And plot it
figure()
plotGrid(Gc, 'facecolor','none')
axis equal tight
title('compositePebiGrid(...)')
drawnow

%% Set pebiGrid Parameters:
% We now use the other wrapper function to create a PEBI-grid using
% distmesh:
eps = 1/12; % This parameter defines the refinement around the wells. The
            % cell size increase exponentially from the wells with 
            %  exp(-(distance from well)/eps);

%% Generate grid
% distmesh will most likely not converge in the maximum number of iterations.
% This is usually not a problem, since the grid most likely is good before
% the convergence requirement is met. 

Gdist = pebiGrid(gS, [1, 1], 'wellLines', well, ...
                'wellGridFactor',wGf/2, 'wellRefinement', true, ...
                'wellEps',eps, 'faultlines', fault,'faultGridFactor',0.8);
Gdist = computeGeometry(Gdist);
%% Plot pebiGrid
figure(); hold on
plotGrid(Gdist,'facecolor','none')
plotGrid(Gdist,Gdist.cells.tag, 'facecolor','b')
centF = Gdist.faces.centroids(Gdist.faces.tag,:);
plot(centF(:,1), centF(:,2),'.','markersize',10)
axis equal tight
title('pebiGrid(...)')



