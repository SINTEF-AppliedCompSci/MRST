%% Impact of Realistic Geologic Models on Simulation of CO2 storage (IGEMS)
% The purpose of the IGEMS project was to investigate how variations in the
% top-surface morphology with a relief amplitude below seismic resolution
% influence CO2 storage capacity. To this end, the project considered fifteen
% different types of top-surface morphologies. For each of these fifteen
% scenarios, except for the flat base case, one hundred different surfaces
% with a 100x100 m resolution were generated stochastically. All data are
% publicly available and can be downloaded from the IGEMS project page
% <http://www.nr.no/IGEMS/>
%
% In this example, we will compute the volumes that potentially are
% available for structural trapping in each of the IGEMS realizations. The
% volumes will be computed using both the cell-based and the node-based
% algorithms. Notice that THIS COMPUTATION WILL BE VERY TIME CONSUMING
% because it altogether involves 1400 realizations that each have 300x600
% cells.

mrstModule add co2lab coarsegrid matlab_bgl;

%% Compute volumes for all realizations
names = dir(fullfile(mrstPath('co2lab'), 'data', 'igems', 'surfaces'));
names = names(3:end);
nrel  = 100;
nvol   = zeros(numel(names),nrel); tn = 0;
cvol   = zeros(numel(names),nrel); tc = 0;
for i=1:numel(names)
   for j=1:nrel
      
      % Load data
      disp(['Realization: ' names(i).name '-' num2str(j)]);
      Gt = topSurfaceGrid(readIGEMSIRAP(names(i).name,j,'save',false));
      
      % Cell-based computation
      t1 = tic;
      res = trapAnalysis(Gt,true);
      cvol(i,j) = sum(res.trap_z(res.traps(res.traps~=0)) - Gt.cells.z(res.traps~=0))*.25e4;
      tc = tc + toc(t1);
      
      % Node-based computation
      t2 = tic;
      res = trapAnalysis(Gt,false);
      nvol(i,j) = sum(res.trap_z(res.traps(res.traps~=0)) - Gt.cells.z(res.traps~=0))*.25e4;
      tn = tn + toc(t2);
   end
end
save igems-results.mat cvol nvol tc tn;

%% Make histograms of volume estimates
% Uses barweb function, download from
% http://www.mathworks.com/matlabcentral/fileexchange/10803-barweb-bargraph-with-error-bars

% Volume plot: cell based
figure
m=mean(cvol,2);  m = reshape([m(1:10); 0; m(11:14)],5,3)'; mc=m([3 1 2],:);
s=std(cvol,0,2); s = reshape([s(1:10); 0; s(11:14)],5,3)'; sc=s([3 1 2],:);
barweb(mc, sc, [],{'flat', 'FMM', 'OSS'}, 'cell based', [], ...
   'trap volume', winter, [], {'No', 'NP1', 'NP2', 'UP1', 'UP2'}, 0, 'plot');
print -depsc2 igems-cbar.eps;

% Volume plot: node based
figure
m=mean(nvol,2);  m = reshape([m(1:10); 0; m(11:14)],5,3)'; mn=m([3 1 2],:);
s=std(nvol,0,2); s = reshape([s(1:10); 0; s(11:14)],5,3)'; sn=s([3 1 2],:);
barweb(mn, sn, [],{'flat', 'FMM', 'OSS'}, 'node based', [], ...
   'trap volume', winter, [], {'No', 'NP1', 'NP2', 'UP1', 'UP2'}, 0, 'plot');
print -depsc2 igems-nbar.eps;

% Volume plot: both in one
figure
m = [mc mn]; m = m(:,[1 6 2 7 3 8 4 9 5 10]);
s = [sc sn]; s = s(:,[1 6 2 7 3 8 4 9 5 10]);
col  = colorcube(8); 
col  = col([1 3:6],:); 
scol = 0.4*col + 0.6*ones(5,3); 
col  = [col; scol];
col  = col([1 6 2 7 3 8 4 9 5 10],:);
barweb(m, s, [], {'flat', 'FMM', 'OSS'}, 'cell/node based', [], ...
   'trap volume', col, [], {'No - cell','No - node','NP1 - cell',...
   'NP1 - node', 'NP2 - cell', 'NP2 - node','UP1 - cell', ...
   'UP1 - node','UP2 - cell','UP2 - node'}, 0, 'plot');
print -depsc2 igems-cnbar.eps;