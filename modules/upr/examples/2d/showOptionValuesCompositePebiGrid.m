%% Example
% In this example we show all optional options for compositePebiGrid2D.

%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2016-2020 Runar Lie Berge. See COPYRIGHT.TXT for details.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}

%% cellConstraints
% This sets the cell constraints in our reservoir. The desired feature of
% these constraints is that the cell centroids of the grid cell trace the
% lines that define the constraints. In practice, we can not guarantee that
% the cell centroids trace the given constraints, instead, the algorithm
% will place the PEBI sites along the constrained lines. However, if the
% grid is fairly "nice", the cell centroids will align with the PEBI sites.
% The cell contstraints are stored as a cell array, each element
% corresponding to one line:
w = {[0.2,0.8;0.5,0.6;0.8,0.8],...
     [0.5,0.2]};
gS = [0.1,0.1];
pdims=[1,1];
G = compositePebiGrid2D(gS, pdims,'cellConstraints',w);

clf
plotGrid(G);
plotGrid(G,G.cells.tag,'facecolor','b')
axis equal tight off
plotLinePath(w,'wo-','linewidth',2,'MarkerFaceColor','w');


%% CCFactor
% The CCFactor sets the relative distance between the sites that trace the
% cellConstraints. If CCFactor=0.5 the distance between the cellConstraint
% sites will be about half the size of the reservoir grid size:
w = {[0.2,0.3;0.8,0.7]};
gS = [0.1,0.1];
pdims=[1,1];
G1 = compositePebiGrid2D(gS, pdims,'cellConstraints',w,'CCFactor',1);
G2 = compositePebiGrid2D(gS, pdims,'cellConstraints',w,'CCFactor',1/2);

figure('Position',[480 340 980 420])
subplot(1,2,1)
plotGrid(G1), title('CCFactor=1'), axis equal tight off
plotGrid(G1,G1.cells.tag,'facecolor','b')
hold on, plotLinePath(w,'wo-','linewidth',2,'MarkerFaceColor','w');
subplot(1,2,2)
plotGrid(G2), title('CCFactor=1/2'), axis equal tight off
plotGrid(G2,G2.cells.tag,'facecolor','b')
hold on, plotLinePath(w,'wo-','linewidth',2,'MarkerFaceColor','w');


%% mlqtMaxLevel
% mlqtMaxLevel sets how many steps of local grid refinement we create
% around each cellConstraint.

w = {[0.2,0.3;0.5,0.5;0.8,0.5]};
gS = [0.1,0.1];
pdims=[1,1];
G1 = compositePebiGrid2D(gS, pdims,'cellConstraints',w,'CCFactor',1/2,'mlqtMaxLevel',1);
G2 = compositePebiGrid2D(gS, pdims,'cellConstraints',w,'CCFactor',1/8,'mlqtMaxLevel',3);

figure('Position',[480 340 980 420])
subplot(1,2,1)
plotGrid(G1); title('mlqtMaxLevel=1'), axis equal tight off
plotGrid(G1,G1.cells.tag,'facecolor','b')
hold on, plotLinePath(w,'wo-','linewidth',2,'MarkerFaceColor','w');
subplot(1,2,2)
plotGrid(G2); title('mlqtMaxLevel=3'), axis equal tight off
plotGrid(G2,G2.cells.tag,'facecolor','b')
hold on, plotLinePath(w,'wo-','linewidth',2,'MarkerFaceColor','w');


%% mlqtLevelSteps
% We can define what distance from a cellConstraint triggers which level of
% mlqt-refinement.
w = {[0.2,0.3;0.5,0.5;0.8,0.5]};
gS = [0.1,0.1];
pdims=[1,1];
lev1 = [0.05,0.1];
lev2 = [0.1,0.2];
G1 = compositePebiGrid2D(gS, pdims,'cellConstraints',w,'CCFactor',1/4, ...
                      'mlqtMaxLevel',2,'mlqtLevelSteps',lev1);
G2 = compositePebiGrid2D(gS, pdims,'cellConstraints',w,'CCFactor',1/4, ...
                      'mlqtMaxLevel',2,'mlqtLevelSteps',lev2);

figure('Position',[480 340 980 420])
subplot(1,2,1)
plotGrid(G1); title('mlqtLevelSteps=[0.05,0.1]'), axis equal tight off
plotGrid(G1,G1.cells.tag,'facecolor','b')
hold on, plotLinePath(w,'wo-','linewidth',2,'MarkerFaceColor','w');
subplot(1,2,2)
plotGrid(G2); title('mlqtLevelSteps=[0.1,0.2]'), axis equal tight off
plotGrid(G2,G2.cells.tag,'facecolor','b')
hold on, plotLinePath(w,'wo-','linewidth',2,'MarkerFaceColor','w');


%% CCRho
% A function that sets the relative distance between the cellConstraint
% sites in the domain.
w = {[0.2,0.3;0.5,0.5;0.8,0.7]};
gS = [0.1,0.1];
pdims=[1,1];
CCRho = @(p) 1 - 0.9*p(:,1);
G1 = compositePebiGrid2D(gS, pdims,'cellConstraints',w);
G2 = compositePebiGrid2D(gS, pdims,'cellConstraints',w,'CCRho',CCRho);

figure('Position',[480 340 980 420])
subplot(1,2,1)
plotGrid(G1); plotGrid(G1,G1.cells.tag,'facecolor','b');
plotGrid(G1,G1.cells.tag,'facecolor','b')
hold on, plotLinePath(w,'wo-','linewidth',2,'MarkerFaceColor','w');
title('CCRho=1'); axis equal tight off
subplot(1,2,2)
plotGrid(G2); plotGrid(G2,G2.cells.tag,'facecolor','b');
plotGrid(G2,G2.cells.tag,'facecolor','b')
hold on, plotLinePath(w,'wo-','linewidth',2,'MarkerFaceColor','w');
title('CCRho=@(p) 1 - 0.9*p(:,1)')
axis equal tight off


%% protLayer
% Adds a protection layer around cellConstraints. The protection sites are placed
% normal along the constraint path
x = linspace(0.2,0.8);
y = 0.5+0.1*sin(pi*x);
w = {[x',y']};
gS = [0.1,0.1];
pdims=[1,1];
G1 = compositePebiGrid2D(gS, pdims, 'cellConstraints', w, ...
                       'interpolateCC', true, 'protLayer', false);
G2 = compositePebiGrid2D(gS, pdims, 'cellConstraints', w, ...
                       'interpolateCC', true, 'protLayer', true);

figure('Position',[480 340 980 420])
subplot(1,2,1)
plotGrid(G1); plotGrid(G1,G1.cells.tag,'facecolor','b')
title('Protection Layer off'), axis equal tight off
subplot(1,2,2)
plotGrid(G2); plotGrid(G2,G2.cells.tag,'facecolor','b')
title('Protection Layer on'), axis equal tight off


%% protD
% This parameter sets the distance from the protection sites to the
% cellConstraint. This is a cellarray of functions, one for each cellConstraint.
w = {[0.2,0.8;0.8,0.8],...
     [0.5,0.2;0.1,0.5]};
gS = [0.08,0.08];
pdims=[1,1];
protD = {@(p)0.01+0.1*p(:,1), @(p) 0.01*ones(size(p,1),1)};
G1 = compositePebiGrid2D(gS, pdims,'cellConstraints',w,'protLayer',true);
G2 = compositePebiGrid2D(gS, pdims,'cellConstraints',w,'protLayer',true,'protD',protD);

figure('Position',[480 340 980 420])
subplot(1,2,1)
plotGrid(G1); plotGrid(G1,G1.cells.tag,'facecolor','b')
title('Default protection distance'), axis equal tight off
subplot(1,2,2)
plotGrid(G2); plotGrid(G2,G2.cells.tag,'facecolor','b')
title('User-prescribed protection distance'), axis equal tight off


%% faceConstraints
% Defines a set of surfaces which should be traced by faces of the grid.
% The faults are stored as a cell array, each element corresponding to one
% fault
f = {[0.2,0.8;0.5,0.65;0.8,0.8],...
     [0.5,0.2;0.1,0.5]};
gS = [0.1,0.1];
pdims=[1,1];
G = compositePebiGrid2D(gS, pdims,'faceConstraints',f);

figure()
plotGrid(G);
plotFaces(G,G.faces.tag,'edgeColor','r','LineWidth',4)
plotLinePath(f,'--ow','linewidth',2,'MarkerFaceColor','w');
axis equal tight off
%% interpolateFC
% The face constraints can either be interpolated or represented exactly.
% the interpolation option should be used if the line path is defined as a
% curve that is discretized by very many line segments. If the line segment
% is interpolated faces of the grid might not represent the face constraints
% exactly at the "kinks" of the curves. If interpolateFC is
% false, each line segment in a face constraint will be represented
% exactly, also where the curve bends.
f = {[0.2,0.8;0.5,0.65;0.8,0.8],...
     [0.5,0.2;0.1,0.5]};
gS = [0.1,0.1];
pdims=[1,1];
G1 = compositePebiGrid2D(gS, pdims,'faceConstraints',f,'interpolateFC',[false; false]);
G2 = compositePebiGrid2D(gS, pdims,'faceConstraints',f,'interpolateFC',[true; true]);

figure('Position',[480 340 980 420])
subplot(1,2,1)
plotGrid(G1); plotFaces(G1,G1.faces.tag,'edgeColor','r','LineWidth',4)
plotLinePath(f,'--ow','linewidth',2,'MarkerFaceColor','w');
title('Fault lines: no interpolation'), axis equal tight off
subplot(1,2,2)
plotGrid(G2); plotFaces(G2,G2.faces.tag,'edgeColor','r','LineWidth',4)
plotLinePath(f,'--ow','linewidth',2,'MarkerFaceColor','w');
title('Fault lines: interpolated'), axis equal tight off


%% FCFactor
% The FCFactor sets the relative distance between the sites that represent
% the constrained surface. If FCFactor=0.5 the sites will
% have spaceing about half the size of the reservoir cells:
f = {[0.2,0.3;0.5,0.5;0.8,0.5]};
gS = [0.1,0.1];
pdims=[1,1];
G1 = compositePebiGrid2D(gS, pdims,'faceConstraints',f,'FCFactor',1);
G2 = compositePebiGrid2D(gS, pdims,'faceConstraints',f,'FCFactor',1/2);

figure('Position',[480 340 980 420])
subplot(1,2,1)
plotGrid(G1); title('FCFactor=1'); axis equal tight off
plotLinePath(f,'--or','linewidth',2,'MarkerFaceColor','w');
subplot(1,2,2)
plotGrid(G2); title('FCFactor=1/2'); axis equal tight off
plotLinePath(f,'--or','linewidth',2,'MarkerFaceColor','w');


%% circleFactor
% The surface sites are generated by setting placing a set of circles along
% the surface path with a distance gS*FCFactor. The radius of the
% circles are gS*FCFactor*circleFactor. The fault sites are place
% at the intersection of these circles.
f = {[0.2,0.3;0.5,0.5;0.8,0.5]};
gS = [0.1,0.1];
pdims=[1,1];
G1 = compositePebiGrid2D(gS, pdims,'faceConstraints',f,'circleFactor',0.55);
G2 = compositePebiGrid2D(gS, pdims,'faceConstraints',f,'circleFactor',0.9);

figure('Position',[480 340 980 420])
subplot(1,2,1)
plotGrid(G1); title('circleFactor=0.55'), axis equal tight off
plotLinePath(f,'--or','linewidth',2,'MarkerFaceColor','w');
subplot(1,2,2)
plotGrid(G2); title('circleFactor=0.9'), axis equal tight off
plotLinePath(f,'--or','linewidth',2,'MarkerFaceColor','w');

%% PolyBdr
% Create a non-square reservoir domain
bdr   = [0,0;-0.5,1;0.5,1.5;1,1];
gs    = [0.1,0.1];
pdims = [1,1]; % can be set to anything
G = compositePebiGrid2D(gs,pdims,'polyBdr',bdr);

figure()
plotGrid(G); title('Polygonal boundary'); axis equal tight off
