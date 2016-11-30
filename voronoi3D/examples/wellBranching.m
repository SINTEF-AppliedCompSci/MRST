%% Create a well branching in 3D
% In this example we will use compositePebiGrid3D to create a grid of a
% reservoir with a well that is branching multiple times.
%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2016 Runar Lie Berge. See COPYRIGHT.TXT for details.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}

%% Set reservoir size
x = 1;
y = 0.5;
z = 0.5;
celldim = [20,10,10];
%% Set well paths
wG = @(p) x/30*ones(size(p,1),1);
wl = {[x/2,y/2,.99*z;x/2,y/2,3*z/4],  ...
      [x/2,y/2,3*z/4; x/4,  y/2, z/2 ; x/6,  y/2,z/6], ...
      [x/2,y/2,3*z/4; 3*x/4,y/2, z/2; 5*x/6,y/2,z/6], ...
      [x/2,y/2,3*z/4;x/2,y/2,z/2], ...
      [x/2,y/2,z/2; 3*x/5,y/2,z/3; 2*x/3,y/2,.01], ...
      [x/2,y/2,z/2; 2*x/5,y/2,z/3;   x/3,y/2,.01]};

%% Create grid
G = compositePebiGrid3D(celldim,[x,y,z], 'wellLines',wl,'wellRho',{wG});
G = computeGeometry(G);
%% plot
color = get(gca,'ColorOrder');
close all
c = G.cells.centroids(:,2) < y/2.01;
plotGrid(G,c)
plotGrid(G,G.cells.tag,'facecolor','red')
view(180,0)
axis equal

