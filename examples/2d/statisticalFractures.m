%% Create a fractrued reservoir
% The following fractures are taken from a dataset in the HFM module in
% MRST. The example full will take at least 8-10 minutes to run, depending
% on your computer. To lower the runtime, we select a subset of the
% fractures by restricting the domain along the y-axis. Likewise, we
% linearize the computation of edge-length distance in distmesh2d.

%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2016 Runar Lie Berge. See COPYRIGHT.TXT for details.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}  
mrstModule add upr

%% Load dataset and possibly select a subset of the fractures
% To select a smaller subset of the fractures, you can set yMax to be
% smaller than the default value of 120.
yMax = 60;
pth = fullfile(mrstPath('upr'), 'datasets', 'statistical_fractures.mat');
load(pth)
l = mat2cell(fl, ones(size(fl,1),1), size(fl,2)).';
l = cellfun(@(c) reshape(c', 2,[])' + 2,l,'un',false);

ix = false(size(l));
for i=1:numel(ix)
  ll=l{i}; 
  ix(i) = max(ll(:,2))<yMax;
end
l = l(ix);

clf;
subplot(2,1,1);
patch([0 35 35 0],[0 0 yMax yMax],ones(4,1),'FaceColor',[.95 .95 1])
plotLinePath(l,'LineWidth',1);
view(90,90), axis tight

%% Set grid parameters.
% We need to refine the grid in some areas, since some fractures are very
% close to each other, but not intersecting. We mark these points and
% introduce a local 
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
faultRho =@(p) min(ones(size(p,1),1), ...
    min(bsxfun(@times, amp,exp(bsxfun(@rdivide, pdist2(p,refPts),eps))),[],2));

% Plot the distance function
G = computeGeometry(cartGrid([70 2*yMax],[35 yMax]));
subplot(2,1,2)
plotCellData(G, faultRho(G.cells.centroids),'EdgeColor','none'); view(90,90)
hold on
plotLinePath(l,'k', 'LineWidth',1);
plot(refPts(:,1),refPts(:,2),'w.','MarkerSize',12);
hold off,
axis tight

%% Generate grid
G = pebiGrid(15,[35,yMax],'faultLines',l,'faultGridFactor',1/50,...
             'circleFactor',0.62,'faultRefinement',true,'faultEps',5,...
             'faultRho', faultRho,'useMrstPebi',true,'linearize', true);

%% Plot the grid
subplot(2,1,1)
plotGrid(G,'FaceColor','none'), axis tight
