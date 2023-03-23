%% Stair-stepped grid
clf
grdecl = readGRDECL(fullfile(getDatasetPath('CaseB4'),'stairstep_36x48.grdecl'));
Gs = processGRDECL(grdecl);
Gs = computeGeometry(Gs);
Kx = logNormLayers(Gs.cartDims, [10 50 400 50 10], ...
   'indices', [1 2 5 8 12 13]); Kx = log10(Kx);

plotCellData(Gs,Kx(Gs.cells.indexMap),'EdgeColor','k');
set(gca,'dataa',[15 20 1])
axis tight off, view(150,30), zoom(1.2)
set(gca,'Clipping','off')

%%
cut_grdecl = cutGrdecl(grdecl, [12 20; 13 23; 1 12]);
k = reshape(Kx,Gs.cartDims);
kx = reshape(k(12:20,13:23,1:12),1,[]);
g = processGRDECL(cut_grdecl);

clf
plotCellData(g,kx(g.cells.indexMap)','EdgeColor','k'),
axis tight off, view(150,50), zoom(1.2)
set(gca,'Clipping','off')

%% Corner-point grid
clf
grdecl = readGRDECL(fullfile(getDatasetPath('CaseB4'),'pillar_36x48.grdecl'));
Gp = processGRDECL(grdecl);
Gp = computeGeometry(Gp);
plotCellData(Gp,Kx(Gp.cells.indexMap),'EdgeColor','k'); 
set(gca,'dataa',[15 20 1])
axis tight off, view(150,30), zoom(1.2)
set(gca,'Clipping','off')

%%
cut_grdecl = cutGrdecl(grdecl, [12 20; 13 23; 1 12]);
g = processGRDECL(cut_grdecl);

clf
plotCellData(g,kx(g.cells.indexMap)','EdgeColor','k'),
axis tight off, view(150,50), zoom(1.2)
set(gca,'Clipping','off')

%% Visualize results
clf
set(gcf,'Position', [300 450 1000 330]);
subplot(1,2,1); 
plotCellData(Gs, log10(Kx(Gs.cells.indexMap)), 'EdgeColor', 'k', 'EdgeAlpha',0.1);
set(gca,'dataaspect',[32 44 2.5])
axis tight off, view(104,28);
% 
subplot(1,2,2);
plotCellData(Gp, log10(Kx(Gp.cells.indexMap)),'EdgeColor', 'k', 'EdgeAlpha',0.1);
set(gca,'dataaspect',[32 44 2.5])
axis tight off, view(104,28);

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
