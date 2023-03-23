%% Create a well branching in 3D
% In this example we will use compositePebiGrid3D to create a grid of a
% reservoir with a well that is branching multiple times.
%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2016 Runar Lie Berge. See COPYRIGHT.TXT for details.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}

%% Set reservoir size
xmax = 1.01;
ymax = 1.01; 
zmax = 1.01;
celldim = [13, 13, 13];

%% Create well sites
wGc = @(p) 1/13 / 4 * ones(size(p, 1), 1);
line1 = [0.5, 0.5, 1; 0.5, 0.5, 3 / 4];
line2 = [0.5, 0.5, 3 / 4; 1 / 5, 0.5, 0.5; 1 / 6, 0.5, 1 / 6];
line3 = [0.5, 0.5, 3 / 4; 4 / 5, 0.5, 0.5; 5 / 6, 0.5, 1 / 6];
line4 = [0.5, 0.5, 3 / 4; 0.5, 0.5, 0.5];
line5 = [0.5, 0.5, 0.5; 3 / 5, 0.5, 1 / 3; 2 / 3, 0.5, 0.1];
line6 = [0.5, 0.5, 0.5; 2 / 5, 0.5, 1 / 3; 1 / 3, 0.5, 0.1];
wl = {line1, line2, line3, line4, line5, line6};

%% Create grid
G = compositePebiGrid3D(celldim, [xmax, ymax, zmax], ...
    'cellConstraints', wl, ...
    'CCRho', {wGc}, ...
    'mlqtMaxLevel', 2, ...
    'mlqtLevelSteps', 0.1);
G = computeGeometry(G);
%% plot
color = get(gca,'ColorOrder');
close all
c = G.cells.centroids(:,2) < 0.5;
plotGrid(G,c)
plotGrid(G,G.cells.tag,'facecolor','red')
view(180,0)
axis equal tight

