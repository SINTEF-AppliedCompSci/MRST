%% Create a fractrued reservoir
% The following fractures are taken from a dataset in the HFM module in MRST
% Before you run this code you probably want to reduce the maximum number
% of iterations in distmesh to around 30-40 or this will take a very long
% time

%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2016 Runar Lie Berge. See COPYRIGHT.TXT for details.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}  

%% load dataset
pth = fullfile(mrstPath('upr'), 'datasets', 'statistical_fractures.mat');
load(pth)
l = mat2cell(fl, ones(size(fl,1),1), size(fl,2));
l = l';
offset = 2;
l = cellfun(@(c) reshape(c', 2,[])' + offset,l,'un',false);

clf;
plotLinePath(l);
%% Set grid parameters.
% We need to refine the grid in some areas since some fractures are very
% close to each other, but not intersecting.
eps = [3,5,1,1,1,2,2,1.5,1,1.5,1.5,1.5,1,1.5,1.5,2.5,1.5];
refPts = [17.4,79.5; 15.0,71.6; 25.9,98.7; 25.8,97.9; ...
          25.7,97.2; 25.7,96.4; 25.5,94.7; 25.3,92.8; ...
          9.1,87.0;  13.0,87.6; 13.1,88.5; 9.0 ,66.2; ...
          28.95,77;  26.0,70.0; 20.8,69.2; 24.9,24.8; 26.8,38.9];
amp = [0.5,0.7,0.3,0.3,0.4,0.4,0.4,.6,.3,.45,.45,.2,.3,0.25,0.55,0.45,0.6];
faultRho =@(p) min(ones(size(p,1),1), ...
                   min(bsxfun(@times, amp,exp(bsxfun(@rdivide, pdist2(p,refPts),eps))),[],2));
%% Generate grid
% Note that this may take some time to execute
G = pebiGrid(15,[35,120],'faultLines',l,'faultGridFactor',1/50,...
             'circleFactor',0.62,'faultRefinement',true,'faultEps',5,...
             'faultRho', faultRho);
           

%% plot grid
figure(); hold on
plotGrid(G)
plotLinePath(l,'--');
axis equal
