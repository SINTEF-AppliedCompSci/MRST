%% Create a fractured reservoir
% The following fractures are taken from a dataset in the HFM module in
% MRST. The example full will take at least 8-10 minutes to run, depending
% on your computer. To lower the runtime, we select a subset of the
% fractures by restricting the domain along the y-axis. Likewise, we
% linearize the computation of edge-length distance in distmesh2d.

%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2016-2020 Runar Lie Berge. See COPYRIGHT.TXT for details.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}

mrstModule add upr

%% Load dataset and possibly select a subset of the fractures
% To select a smaller subset of the fractures, you can set yMax to be
% smaller than the default value of 120.
yMax = 120;
pth = fullfile(mrstPath('upr'), 'datasets', 'statistical_fractures.mat');
load(pth)
lines = mat2cell(fl, ones(size(fl,1),1), size(fl,2)).';
lines = cellfun(@(c) reshape(c', 2,[])' + 2,lines,'un',false);

ix = false(size(lines));
for i=1:numel(ix)
  ll=lines{i}; 
  ix(i) = max(ll(:,2))<yMax;
end
lines = lines(ix);

clf;
subplot(2,1,1);
patch([0 35 35 0],[0 0 yMax yMax],ones(4,1),'FaceColor',[.95 .95 1])
plotLinePath(lines,'LineWidth',1);
view(90,90), axis tight

%% Set grid parameters.
% We need to refine the grid in some areas, since some fractures are very
% close to each other, but not intersecting. We mark these points and
% introduce local grid refinement by setting up a grid density function
% whose values decay exponentially as we approach any of the points of
% interest.
eps = [3,5,1,1,1,2,2,1.5,1,1.5,1.5,1.5,1,1.5,1.5,2.5,1.5];
refPts = [17.4,79.5; 15.0,71.6; 25.9,98.7; 25.8,97.9; ...
          25.7,97.2; 25.7,96.4; 25.5,94.7; 25.3,92.8; ...
          9.1,87.0;  13.0,87.6; 13.1,88.5; 9.0 ,66.2; ...
          28.95,77;  26.0,70.0; 20.8,69.2; 24.9,24.8; 26.8,38.9];
amp = [0.5,0.7,0.3,0.3,0.4,0.4,0.4,.6,.3,.45,.45,.2,.3,0.25,0.55,0.45,0.6];
ix = false(size(eps));
for i=1:numel(ix) 
  ix(i) = refPts(i,2)<yMax;
end
[eps,amp,refPts] = deal(eps(ix),amp(ix),refPts(ix,:));
FCRho =@(p) min(ones(size(p,1),1), ...
    min(bsxfun(@times, amp,exp(bsxfun(@rdivide, pdist2(p,refPts),eps))),[],2));

% Plot the distance function
G = computeGeometry(cartGrid([70 2*yMax],[35 yMax]));
subplot(2,1,2)
plotCellData(G, FCRho(G.cells.centroids),'EdgeColor','none'); view(90,90)
hold on
plotLinePath(lines,'k', 'LineWidth',1);
plot(refPts(:,1),refPts(:,2),'r.','MarkerSize',12);
patch([0 35 35 0],[0 0 yMax yMax],ones(4,1),'FaceColor','none')
colormap([parula.^2; 1 1 1])
hold off,
axis tight

%% Generate grid
% We generate a grid where the fractures are traced by faces of the grid.
% This will take a few minutes since the number of constraints is quite
% large.
G = pebiGrid2D(15,[35,yMax],'faceConstraints',lines,'FCFactor',1/50,...
             'circleFactor',0.55,'FCRefinement',true,'FCEps',5,...
             'FCRho', FCRho,'useMrstPebi',true,'linearize', true);

%% Plot the grid
clf
set(gcf,'Position',[500 370 878 408]);
subplot(3,1,1:2)
patch([0 35 35 0],[0 0 yMax yMax],ones(4,1),'FaceColor',[.95 .95 1]), 
plotGrid(G,'FaceColor','none','edgecolor',[.6 .6 .6])
h=plotLinePath(lines,'LineWidth',1.5);
mrstModule add book
col = tatarizeMap(numel(h));
for i=1:numel(h), set(h(i),'Color',col(i,:)); end
hold on
patch([23.5 25.5 25.5 23.5],[17 17 28 28],ones(4,1), ...
    'LineWidth',1,'EdgeColor','b','FaceColor','none')
patch([24.6 27 27 24.6],[65 65 74 74],ones(4,1), ...
    'LineWidth',1,'EdgeColor','r','FaceColor','none')
hold off
view(90,90), axis tight off
ax1 = gca;

% Zoom 1
subplot(3,2,5)
copyobj(get(ax1,'Children'),gca);
axis tight off, axis([23.5 25.5 17 28]); view(90,90)
ax=gca; ax.Position = ax.Position + [0 0 .03 .04];

% Zoom 2
subplot(3,2,6)
copyobj(get(ax1,'Children'),gca);
axis tight off, axis([24.6 27 65 74]); view(90,90)
ax=gca; ax.Position = ax.Position + [-.03 0 .03 .04];

