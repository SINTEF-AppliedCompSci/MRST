%% Exmple
% This example contains three faults that intersects. The example is
% inspired by Ding et Al.:
% X. Y. Ding and L. S. K. Fung. An unstructured Gridding Method for
% Simulating Faulted Reservoirs polulated with Complex Wells. In
% Proceedings of the SPE Reservoir Simulation Symposium, Houston, Texas,
% USA, February 2015. doi: 10.2118/173243-MS
%
%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2016 Runar Lie Berge. See COPYRIGHT.TXT for details.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}  
%% Set fault
l = {[0.1,0.42; 0.4,.55; 0.7,0.65], ...
     [0.8,0.13; 0.6,0.4; 0.55,0.6],...
     [0.42,1.08; 0.45,0.9; 0.5,0.8; 0.58,0.6]};
   
gS = [1/30,1/30];
   
%% Create grid
% We use the routine compositePebiGrid2D to create the grid. We use the
% preset options. 
G = compositePebiGrid2D(gS,[1,1.15],'faceConstraints',l, 'useMrstPebi', true);

%% Plot grid
figure(1); clf; hold on
plotGrid(G,'facecolor','none');
plotLinePath(l,'-','color','k','linewidth',1);
axis equal

%% Plot zoomed
% We plot a zoomed view of the intersection
figure(2); clf; hold on
plotGrid(G,'facecolor','none');
plotLinePath(l,'-','color','k','linewidth',1.5);
axis equal off tight
axis([0.4,0.7,0.4,0.7 + (0.7-0.4)*.15])

