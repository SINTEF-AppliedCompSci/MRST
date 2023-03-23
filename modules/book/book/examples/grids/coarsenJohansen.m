%% Coarsen Real Models: the Johansen Formation
% In this script we will look at how to coarsen real model: a
% sector model of the Johansen aquifer which was considered as a possible
% candidate for large-scale geological storage of CO2. The data can be
% downloaded from
% http://www.sintef.no/Projectweb/MatMorA/Downloads/Johansen/

mrstModule add coarsegrid;

dpath = getDatasetPath('johansen');
sector = fullfile(dpath, 'NPD5');
filename = [sector, '.grdecl'];
G = processGRDECL(readGRDECL(filename));
G = computeGeometry(G);
K = reshape(load([sector, '_Permeability.txt'])', prod(G.cartDims), []);
K = K(G.cells.indexMap);

%%
% First, we make a partition with a coarsening factor four in each lateral
% direction. In the vertical direction, the model consists of three
% different formations that can be distinguised by the permeability values:
% the Dunlin shale has K<=0.01mD, the Amundsen shale has K<0.1mD, whereas
% the Johansen sandstone has K>1mD. For illustration purposes, we only keep
% one block in the vertical direction for each of the formations.
pK = 2*ones(size(K));  % Johansen
pK(K<=0.1) = 3;        % Amundsen
pK(K<=0.01)= 1;        % Dunlin
pC = partitionUI(G, [G.cartDims(1:2)/4 1]);
[~,~,p] = unique([pK, pC], 'rows');
p = processPartition(G,p);

clf
plotCellData(G,log10(K),'EdgeColor','k','EdgeAlpha',.4); view(3)
outlineCoarseGrid(G,p,'FaceColor','none','EdgeColor','k','LineWidth',1.5);
axis tight off
colormap((2*jet(16)+ones(16,3))/3);

%%
% Make a new coarse grid in which we keep the vertical resolution, which is
% typically what one would do if the model is to be used to simulate CO2
% storage. The resulting coarse model will have many coarse blocks that
% have more than six connections. We loop through the blocks in the top
% layer to show how they look like
pK = 2*ones(size(K)); pK(K<=0.1) = 3; pK(K<=0.01)= 1;
pC = partitionUI(G, G.cartDims./[4 4 1]);
[~,~,p] = unique([pK, pC], 'rows');
p = processPartition(G,p);
CG = generateCoarseGrid(G,p);
cn = diff(CG.cells.facePos);

figure('Position',[0 60 560 820])
subplot(2,1,1);
plotGrid(G,'FaceColor','none','EdgeAlpha',.1);
plotGrid(CG,cn>6);

h1=[];
val = cumsum(ones(G.cartDims),3); val=val(G.cells.indexMap);
kn = zeros(size(cn)); kn(p) = val;
ind = find(cn>6 & kn==1);
for i=1:numel(ind)
   subplot(2,1,1); delete(h1);
   h1 = plotGrid(CG,ind(i),'FaceColor','r');
   
   subplot(2,1,2); cla;
   plotGrid(G,p==ind(i)); 
   plotGrid(CG,ind(i),'FaceColor','none','LineWidth',2); 
   title(['Block number: ' num2str(ind(i))]);
   view(3); drawnow;
   
end

%%
block = [89 139 143 195 286 466];
ind = [1 5:8 4];
c   = colorcube(9);
col = (2*c(1:6,:)+ones(6,3))/3;
figure('position',[0 60 1000 400]);
subplot(2,4,2:3)
plotGrid(G,'FaceColor','none','EdgeAlpha',.1); 
axis off; zoom(1.3);
for i=1:6
   subplot(2,4,2:3);
   plotGrid(CG,block(i),'FaceColor',col(i,:));
   subplot(2,4,ind(i));
   plotGrid(G,p==block(i),'FaceColor',col(i,:));
   plotGrid(CG,block(i),'FaceColor','none','LineWidth',2);
   view(3); set(gca,'Clipping','off'), axis tight off 
   zoom(1.3)
end

%{
Copyright 2009-2018 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}
