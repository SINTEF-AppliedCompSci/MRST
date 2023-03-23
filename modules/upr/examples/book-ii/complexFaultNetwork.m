%% Example
% This example is inspired by a reservoir a paper by Branets et Al.:
% L. Branets, S.S. Ghai, S. L. Lyons, and X.-H. Wu. Efficient and Accurate
% Reservoir Modeling Using Adaptive Gridding with Global Scale Up. In
% Proceedings of the SPE Reservoir Simulation Symposium, The Woodlands,
% Texas, Jan. 2009. doi: 10.2118/118946-MS.

%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2016-2020 Runar Lie Berge. See COPYRIGHT.TXT for details.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}  

%% Load geometry.
% The faults was drawn in inkscape, and then we extracted the coordinates
% from the .svg file. The dataset contains two variables; fault which is a
% cell array of all faults. The second variable is bdr which is the 
% boundary of the reservoir. bdr is a set of coordinates for the vertices
% in a polygon. The vertices must be ordered counterclockwise or clockwise.
% We can not use the wrapping function pebiGrid2D since the boundary is not a
% square, and therefore have to create the grid "manually". 
pth = fullfile(mrstPath('upr'), 'datasets', 'gridBranets.mat');
load(pth);
figure('Position',[480 340 900 420])
subplot(1,3,1)
col = get(gca, 'colororder');
hold on
plot([bdr(:,1);bdr(1,1)],[bdr(:,2);bdr(1,2)],'k');
plotLinePath(fault, 'color', col(2,:));
axis equal off tight;

%% calculate fault intersections
[faultSplit, fCut, ~] = splitAtInt2D(fault, {});

%% Create fault sites
fGs = max(bdr(:))/50;
F = surfaceSites2D(faultSplit, fGs,'fCut',fCut,'interpolateFC',true);
% Remove tip sites that violate the fault condition
F.t.pts = surfaceSufCond2D(F.t.pts, F);
% Remove tip sites outside domain
innside = inpolygon(F.t.pts(:, 1), F.t.pts(:, 2), bdr(:, 1), bdr(:, 2));
F.t.pts = F.t.pts(innside, :);

%% Create reservoir sites
% We create a set of random reservoir sites. Then we remove any outside our
% domain.
n = 1500;
pInit = bsxfun(@plus,bsxfun(@times, rand(n,2),(max(bdr) - min(bdr))), min(bdr));
keep = inpolygon(pInit(:,1),pInit(:,2),bdr(:,1), bdr(:,2));
pInit = pInit(keep,:);

%% Plot initial sites
pInit = removeConflictPoints(pInit, F.f.pts, F.f.Gs);
subplot(1, 3, 1);
plot(pInit(:,1),pInit(:,2), '.k','MarkerSize', 5);
plot(F.f.pts(:, 1), F.f.pts(:, 2), '.', 'color', col(1,:),'markersize',5);
plot(F.t.pts(:, 1), F.t.pts(:, 2), '.', 'color', col(3,:), 'markersize',5);

%% Create CVD.
% We now create the CVD by minimizing the CVD energy function using the
% L-BFGS agorithm. We do this by calling the wrapping function CPG2D.
% This takes some time (approximately 90 sec), so we have saved
% the output in a file: 'datasets/reservoirWithComplexFaultNetwork.mat'.
% In this file maxIt was set to 10. The CVD algorithm was not fully
% converged, but still gives a quite nice mesh.
% If you wish to run the optimization yourself, uncomment the following
% code:
%[G, sites] = CPG2D(pInit, bdr,'fixedPts', [F.f.pts; F.t.pts],'maxIt', 10);
pth = fullfile(mrstPath('upr'), 'datasets', 'reservoirWithComplexFaultNetwork.mat');
load(pth)

%% Plot grid
subplot(1,3,2), hold on
plotGrid(G,'facecolor','none');
plotLinePath(faultSplit, 'color','k','linewidth',1);
axis equal off tight

%% Create grid using Delaunay optimization
% An alternative to use the CVD optimization is to use the Delaunay
% optimization. This can be used by calling the wrapper function pebiGrid2D

[Gd,ss,F] = pebiGrid2D(fGs, [], 'faceConstraints', fault, ...
                       'polyBdr', bdr, ...
                       'interpolateFC', true);
subplot(1, 3, 3);
plotGrid(Gd,'facecolor','none');
plotLinePath(fault, 'color','k','linewidth',1);
axis equal off tight
