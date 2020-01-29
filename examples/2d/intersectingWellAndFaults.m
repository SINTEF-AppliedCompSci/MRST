%% Example
% This script contains an example of a single curved well path intersected
% by several straight faults. We use the two wrapper-functions
% compositePebiGrid(..) and pebiGrid(..) to create PEBI grids conforming
% to the faults and the well path.

%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2016-2020 Runar Lie Berge. See COPYRIGHT.TXT for details.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}  

%% Set well and fault paths
% Store well and fault paths two cell arrays, in which each cell represents
% a single fracture or well paths. These are in turn given as arrays of
% of points that describe continuous line segments.

% We start by creating a well
x     = linspace(0.1, 0.8, 10);
y     = 0.2 + x.^2;
well  = {[x' , y']};

% And four faults crossing the well.
fault   = {[0.3,0.1;0.1,0.3],...        
          [0.3,0.1;0.4,0.7],...
          [0.6,0.2;0.3,0.7],...
          [0.8,0.4;0.5,0.8]};

% We now plot the well and fault to see how they look like
clf,set(gcf,'Position',[480 340 980 420]);
subplot(1,2,1), ax=gca;
patch([0 1 1 0],[0 0 1 1],ones(4,1),'FaceColor',[.95 .95 1]);
hold all
plotLinePath(well,'-o', 'LineWidth',2,'MarkerSize',3);
plotLinePath(fault,'-o','LineWidth',2,'MarkerSize',3);
title('Well path and fractures')
axis equal tight, axis ([0,1,0,1])
%subplot(1,2,2), copyobj(get(ax,'Children'),gca);

%% Construct Voronoi grid using compositePebiGrid
% Before we call the gridding fuctions we set some parameters.
gS  = 1/24; % The grid size
wGf = 0.5;  % The relative size of well cells compared to gS
fGf = 0.5;  % The relative size of fault cells compared to gS
nRs = 1;    % Number of refinement steps towards the well

% We can now create the composite Pebi grid:
Gc = compositePebiGrid([gS,gS], [1, 1], ...
                       'wellLines', well, 'wellGridFactor',wGf, ...
                       'faultLines',fault, 'faultGridFactor', fGf,...
                       'mlqtMaxLevel', nRs);
% Plot composite grid
subplot(1,2,1)
plotGrid(Gc, 'facecolor','none')
axis equal tight off
title('compositePebiGrid(...)')
drawnow

%% Construct Voronoi grid using pebiGrid
% We now use the other wrapper function to create a PEBI-grid using
% distmesh:
eps = 1/12; % This parameter defines the refinement around the wells. The
            % cell size increase exponentially from the wells with 
            %  exp(-(distance from well)/eps);

% distmesh will most likely not converge in the maximum number of
% iterations. This is usually not a problem, since the grid most likely is
% good before the convergence requirement is met.
Gdist = pebiGrid(gS, [1, 1], 'wellLines', well, ...
                'wellGridFactor',wGf/2, 'wellRefinement', true, ...
                'wellEps',eps, 'faultlines', fault,'faultGridFactor',0.8);
Gdist = computeGeometry(Gdist);

% Plot pebiGrid
subplot(1,2,2), hold on
plotGrid(Gdist,'facecolor',[.95 .95 1])
plotGrid(Gdist,Gdist.cells.tag, 'facecolor','b')
centF = Gdist.faces.centroids(Gdist.faces.tag,:);
plot(centF(:,1), centF(:,2),'.k','markersize',8)
axis equal tight off
title('pebiGrid(...)')



