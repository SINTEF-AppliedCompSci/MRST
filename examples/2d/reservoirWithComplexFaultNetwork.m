%% Example
% This example is inspired by a reservoir a paper by Branets et Al.:
% L. Branets, S.S. Ghai, S. L. Lyons, and X.-H. Wu. Efficient and Accurate
% Reservoir Modeling Using Adaptive Gridding with Global Scale Up. In
% Proceedings of the SPE Reservoir Simulation Symposium, The Woodlands,
% Texas, Jan. 2009. doi: 10.2118/118946-MS.

%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2016 Runar Lie Berge. See COPYRIGHT.TXT for details.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}  

%% load Geometry.
% The faults was drawn in inkscape, and then we extracted the coordinates
% from the .svg file. The dataset contains two variables; fault which is a
% cell array of all faults. The second variable is bdr which is the 
% boundary of the reservoir. bdr is a set of coordinates for the vertices
% in a polygon. The vertices must be ordered counterclockwise or clockwise.
% We can not use the wrapping function pebiGrid since the boundary is not a
% square, and therefore have to create the grid "manually". 
pth = fullfile(mrstPath('upr'), 'datasets', 'gridBranets.mat');
load(pth);
color = get(gca,'ColorOrder');
figure(); hold on
plot([bdr(:,1);bdr(1,1)],[bdr(:,2);bdr(1,2)],'k');
plotLinePath(fault,'color',color(2,:));
axis equal;
%% calculate fault intersections
[fault, fCut, ~] = splitAtInt(fault, {});

%% Create fault sites
fGs = max(max(bdr))/70;
F = createFaultGridPoints(fault, fGs,'fCut',fCut);
plot(F.f.pts(:,1), F.f.pts(:,2),'.','markersize',15);

%% Create reservoir sites
% We crate a set of random reservoir sites. Then we remove any outside our
% domain.
n = 1500;
pInit = bsxfun(@plus,bsxfun(@times, rand(n,2),(max(bdr) - min(bdr))), min(bdr));
keep = inpolygon(pInit(:,1),pInit(:,2),bdr(:,1), bdr(:,2));
pInit = pInit(keep,:);

pInit = removeConflictPoints2(pInit,F.f.pts,F.f.Gs);
%% Create CVD.
% We now create the CVD by minimizing the CVD energy function using the
% L-BFGS agorithm. We do this by calling the wrapping function createCVD.
% This takes some time, so we have saved
% the output in a file: 'datasets/reservoirWithComplexFaultNetwork.mat'.
% In this file maxIt was set to 200, but you get good results even after 10
% iterations.
% If you wish to run the optimization yourself, uncomment the following
% code:
%G = CVD2D(pInit, bdr,'fixedPts', F.f.pts,'maxIt',10);
pth = fullfile(mrstPath('upr'), 'datasets', 'reservoirWithComplexFaultNetwork.mat');
load(pth)

%% Plot grid
figure(); hold on
plotGrid(G,'facecolor','none');
plotLinePath(fault,'color','k','linewidth',1);
axis equal off tight





